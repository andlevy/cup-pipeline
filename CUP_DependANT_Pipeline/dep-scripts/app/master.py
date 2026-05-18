from app.patient_processing import Patient
from app.tools.logger import create_logger
import json

import argparse

logger = create_logger()


def predict_essential_drugs(args, patient_id):
    """The args are mutations_file,expressions_file,disease,keep. Predict uses the mutations_file and expressions_file
    to generate a list of vulnerabilities. Will be modified to include the drugs
    associated with them. This is very much an actively changing function."""

    mutations_file, expressions_file, disease, location = args

    logger.info(f'''mutations_file is {mutations_file},
                    expressions_file is {expressions_file}, 
                    disease is {disease}, 
                    location is {location}''')

    patient = Patient(mutations_file, expressions_file, disease, location, patient_id)

    print('0% - Acquiring mutations...')
    mutations = patient.read_mutations_from_file()

    print('5% - Acquiring expressions...')
    expressions = patient.read_expressions_from_file()

    print('10% - Identifying lofs and gofs...')
    lofs, gofs, gof_rep = patient.identify_lofs_and_gofs(mutations)

    print('25% - Calculating network features...')
    network_features = patient.calculate_network_features(expressions, lofs, gofs)

    print('40% - Predicting essential genes...')
    predicted_essential_genes = patient.predict_essential_genes(network_features)

    print('70% - Finding associated drugs...')
    essential_drugs_df = patient.find_associated_drugs(predicted_essential_genes)

    print('80% - Searching clinical trials...')
    clinical_trials = patient.get_clinical_trials(essential_drugs_df)
    with open(patient.paths.clinical_trials, 'w') as f:
        json.dump(clinical_trials, f, indent=4)

    print('95% - Summarising information...')
    summary = patient.get_summary(list(gof_rep), list(lofs), essential_drugs_df, clinical_trials)

    return summary, clinical_trials


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict essential drugs.")
    parser.add_argument('--mutations', type=str, help='Mutations File')
    parser.add_argument('--expressions', type=str, help='Expressions file')
    parser.add_argument('--disease', type=str, help='Cancer type')
    parser.add_argument('--location', type=str, help='Patient location')
    parser.add_argument('--patient_id', type=str, help='A unique identifier for that patient')

    parsed_args = parser.parse_args()

    input_args = (parsed_args.mutations, parsed_args.expressions, parsed_args.disease, parsed_args.location)
    predict_essential_drugs(input_args, parsed_args.patient_id)
