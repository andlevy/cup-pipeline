import threading
from pathlib import Path
from flask import jsonify, request, url_for
from dummy_app import app
import time
import json
import random
from multiprocessing import Manager

manager = Manager()
tasks = manager.dict()


class Celery:

    def __init__(self):
        self.state = "PENDING"
        self.meta = {"status": "Setting up task..."}
        self.id = random.randint(0, 10000)
        self.task_response = None

    def update_state(self, state, meta):
        self.state = state
        self.meta = meta

    def finish_task(self, task_response):
        if self.task_response:
            raise NotImplementedError("Can't change task after finished")
        self.task_response = task_response
        self.update_state("SUCCESS", {"status": "Task successfully completed"})

    def _execute_task_in_thread(self, method):
        """
        This method is run in a separate thread.
        It executes the user-defined 'method' (e.g., 'run') and handles
        updating the task's state and finishing it.
        """
        try:
            # The 'method' (e.g., 'run') needs to accept the Celery instance
            # so it can call update_state on it.
            result = method(self)
            self.finish_task(result)
        except Exception as e:
            # If an error occurs during task execution, update the state to FAILURE
            self.update_state("FAILURE", {"status": f"Task failed: {str(e)}"})
            print(f"[Task {self.id}] An error occurred: {e}")

    @staticmethod
    def delay(method):
        """
        Launches a task asynchronously in a new thread.
        Returns the Celery task object immediately, allowing the caller
        to poll its status later.
        """
        task = Celery()
        tasks[task.id] = task
        # Create and start a new thread to execute the task
        thread = threading.Thread(target=task._execute_task_in_thread, args=(method,))
        thread.daemon = True  # Allow the main program to exit even if this thread is still running
        thread.start()
        print(f"Launched task {task.id} asynchronously.")
        return task


def run(self):

    self.update_state(state='PROGRESS', meta={'current': 0, 'total': 100, 'status': 'Acquiring mutations...'})
    time.sleep(10)
    self.update_state(state='PROGRESS', meta={'current': 5, 'total': 100, 'status': 'Acquiring expressions...'})
    time.sleep(10)
    self.update_state(state='PROGRESS', meta={'current': 10, 'total': 100, 'status': 'Identifying lofs and gofs...'})
    time.sleep(10)
    self.update_state(state='PROGRESS',
                      meta={'current': 25, 'total': 100, 'status': 'Calculating network features...'})
    time.sleep(10)
    self.update_state(state='PROGRESS',
                      meta={'current': 40, 'total': 100, 'status': 'Predicting essential genes...'})
    time.sleep(10)
    self.update_state(state='PROGRESS', meta={'current': 70, 'total': 100, 'status': 'Finding associated drugs...'})
    time.sleep(10)
    self.update_state(state='PROGRESS',
                      meta={'current': 80, 'total': 100, 'status': 'Searching clinical trials...'})
    time.sleep(10)
    self.update_state(state='PROGRESS', meta={'current': 95, 'total': 100, 'status': 'Summarising information...'})

    with open(Path(__file__).parent / "summary.json", "r+") as f:
        summary = json.load(f)
    with open(Path(__file__).parent / "clinical_trials.json", "r+") as f:
        clinical_trials = json.load(f)
    return {"summary": summary,
            "clinical_trials": clinical_trials}


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

    task = Celery.delay(run)

    return jsonify({
        "message": "Task initiated successfully",
        "task_id": task.id,
        "status_url": url_for('get_task_status', task_id=task.id, _external=True)
    }), 202


@app.route('/tasks/<string:task_id>/status', methods=['GET'])
def get_task_status(task_id):
    task = tasks[task_id]
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
        error_info = task.info
        error_type = error_info.get('error_type', 'UnknownError')
        error_message = error_info.get('message', 'An unknown error occurred during processing.')
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
    task = tasks[task_id]
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
    task = tasks[task_id]
    if task.state == 'SUCCESS':
        clinical_trials = task.get().get("clinical_trials")
        if not clinical_trials:
            response = {"status": "There are no clinical trials",
                        "clinical_trials": [],
                        "state": task.state}
        elif gene:
            if gene not in clinical_trials.keys():
                response = {"status": "There are no clinical trials for that gene",
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
                    response = {"clinical_trials": gene_trials, "state": task.state}
        else:
            response = {"clinical_trials": clinical_trials, "state": task.state}
        return jsonify(response), 200
    elif task.state == 'FAILURE':
        return jsonify({"status": "Task failed", "details": str(task.info), "state": task.state}), 500
    else:
        return jsonify({"status": "Task not yet complete or pending", "state": task.state}), 202
