import json

from flask import request, jsonify, url_for

from app import app, celery
from app.patient_processing import Patient
from app.tools.exceptions import BackendError, FrontendError, IncorrectInputError
from app.tools.logger import create_logger

logger = create_logger()


@celery.task(bind=True)
def run(self,
        mutations_file,
        expressions_file,
        disease,
        location):
    task_id = self.request.id
    try:
        patient = Patient(mutations_file, expressions_file, disease, location, task_id)

        self.update_state(state='PROGRESS', meta={'current': 0, 'total': 100, 'status': 'Acquiring mutations...'})
        mutations = patient.read_mutations_from_file()

        self.update_state(state='PROGRESS', meta={'current': 5, 'total': 100, 'status': 'Acquiring expressions...'})
        expressions = patient.read_expressions_from_file()

        self.update_state(state='PROGRESS',
                          meta={'current': 10, 'total': 100, 'status': 'Identifying lofs and gofs...'})
        lofs, gofs, gof_rep = patient.identify_lofs_and_gofs(mutations)

        self.update_state(state='PROGRESS',
                          meta={'current': 25, 'total': 100, 'status': 'Calculating network features...'})
        network_features = patient.calculate_network_features(expressions, lofs, gofs)

        self.update_state(state='PROGRESS',
                          meta={'current': 40, 'total': 100, 'status': 'Predicting essential genes...'})
        predicted_essential_genes = patient.predict_essential_genes(network_features)

        self.update_state(state='PROGRESS', meta={'current': 70, 'total': 100, 'status': 'Finding associated drugs...'})
        essential_drugs_df = patient.find_associated_drugs(predicted_essential_genes)

        self.update_state(state='PROGRESS',
                          meta={'current': 80, 'total': 100, 'status': 'Searching clinical trials...'})
        clinical_trials = patient.get_clinical_trials(essential_drugs_df)
        with open(patient.paths.clinical_trials, 'w') as f:
            json.dump(clinical_trials, f, indent=4)

        self.update_state(state='PROGRESS', meta={'current': 95, 'total': 100, 'status': 'Summarising information...'})
        summary = patient.get_summary(list(gof_rep), list(lofs), essential_drugs_df, clinical_trials)

        return {"summary": summary.get_results(),
                "clinical_trials": clinical_trials}
    except BackendError:
        logger.exception(f"Task {task_id} failed with BackendError")
        raise
    except FrontendError:
        logger.exception(f"Task {task_id} failed with FrontendError")
        raise
    except IncorrectInputError:
        logger.exception(f"Task {task_id} failed with IncorrectInputError")
        raise
    except Exception as e:
        logger.exception(f"Task {task_id} failed with an unexpected exception")
        raise BackendError from e


@app.route('/')
def index():
    return jsonify({
        "message": "Website backend active"
    }), 200


@app.route('/process_data', methods=['POST'])
def process_data():
    data = request.get_json()

    if not data:
        return jsonify({"error": "Request body must be JSON"}), 400

    required_params = ['mutations_file_name', 'expressions_file_name']
    if not all(param in data for param in required_params):
        return jsonify({"error": "Missing one or more required parameters"}), 400

    mutations_file_name = data.get('mutations_file_name')
    expressions_file_name = data.get('expressions_file_name')
    disease = data.get('disease') if 'disease' in data else ""
    location = data.get('location') if 'location' in data else ""

    try:
        task = run.delay(mutations_file_name, expressions_file_name, disease, location)
    except Exception as e:
        logger.exception("Failed to initialise task")
        return jsonify({
            "error": "Failed to initialise celery task, reason unknown"
        }), 500

    return jsonify({
        "message": "Task initiated successfully",
        "task_id": task.id,
        "status_url": url_for('get_task_status', task_id=task.id, _external=True)
    }), 202


@app.route('/tasks/<string:task_id>/status', methods=['GET'])
def get_task_status(task_id):
    task = celery.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'status': 'Task is pending or not found...'
        }
    elif task.state == 'PROGRESS':
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 100),
            'status': task.info.get('status', '')
        }
    elif task.state == 'SUCCESS':
        response = {
            'state': task.state,
            'status': 'Task has completed successfully',
            'result_available': True
        }
    elif task.state == "FAILURE":
        exc = task.result
        error_type = type(exc).__name__
        error_message = str(exc)
        response = {
            'state': task.state,
            'error_type': error_type,
            'status': error_message
        }
        # Map custom error types to HTTP status codes for the frontend
        if error_type == "IncorrectInputError":
            return jsonify(response), 400
        elif error_type in ["BackendError", "FrontendError"]:
            return jsonify(response), 500
        else:  # Fallback for truly unknown error types
            return jsonify(response), 500
    else:
        response = {
            'state': task.state,
            'status': 'Unknown state'
        }
    return jsonify(response)


@app.route('/tasks/<string:task_id>/summary', methods=['GET'])
def get_summary(task_id):
    task = celery.AsyncResult(task_id)
    if task.state == 'SUCCESS':
        return jsonify({"status": "Task succeeded",
                        "state": task.state,
                        "summary": task.get().get("summary")})
    elif task.state == 'FAILURE':
        return jsonify({"status": "Task failed",
                        "details": str(task.info),
                        "state": task.state}), 500
    else:
        return jsonify({"status": "Task not yet complete or pending",
                        "state": task.state}), 202


@app.route('/tasks/<string:task_id>/clinical_trials/', defaults={'gene': None, 'drug': None}, methods=['GET'])
@app.route('/tasks/<string:task_id>/clinical_trials/<string:gene>', defaults={'drug': None}, methods=['GET'])
@app.route('/tasks/<string:task_id>/clinical_trials/<string:gene>/<string:drug>', methods=['GET'])
def get_clinical_trials(task_id, gene, drug):
    task = celery.AsyncResult(task_id)
    if task.state == 'SUCCESS':
        clinical_trials = task.get().get("clinical_trials")
        if not clinical_trials:
            response = {"status": "There are no clinical trials",
                        "clinical_trials": [],
                        "state": task.state}
        elif gene:
            if gene not in clinical_trials.keys():
                response = {"status": f"There are no clinical trials for {gene}",
                            "clinical_trials": [], "state": task.state}
            else:
                gene_trials = clinical_trials.get(gene)
                if drug:
                    if drug not in gene_trials.keys():
                        response = {"status": f"No clinical trials for {gene}, {drug} combination",
                                    "clinical_trials": [], "state": task.state}
                    else:
                        response = {"status": f"Clinical trials for {gene}, {drug}",
                                    "clinical_trials": gene_trials.get(drug), "state": task.state}
                else:
                    response = {"clinical_trials": gene_trials,
                                "state": task.state,
                                "status": f"Clinical trials for {gene}"}
        else:
            response = {"clinical_trials": clinical_trials,
                        "state": task.state,
                        "status": "All clinical trials"}
        return jsonify(response), 200
    elif task.state == 'FAILURE':
        return jsonify({"status": "Task failed", "details": str(task.info), "state": task.state}), 500
    else:
        return jsonify({"status": "Task not yet complete or pending", "state": task.state}), 202
