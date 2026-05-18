import shutil

import pandas as pd

from app.calculate_network_features import CalculateNetworkFeatures
from app.clinical_trials.get_clinical_trials import GetClinicalTrials
from app.expressions_parser import Expressions
from app.genes_to_drugs import GenesToDrugs
from app.lofs_and_gofs.lofs_and_gofs import IdentifyMutations
from app.mutations_parser import Mutations
from app.predict_essential_genes import Predict
from app.summary import Summary
from app.tools.file_paths import Paths
from app.tools.logger import create_logger
from app.tools.exceptions import BackendError, FrontendError

logger = create_logger()


class Patient:

    def __init__(self,
                 mutations_file: str,
                 expressions_file: str,
                 disease: str,
                 location: str,
                 task_id: str):
        self.paths: Paths = Paths(task_id,
                                  mutations_file.split('/')[-1],
                                  expressions_file.split('/')[-1])

        # move files to a unique task_id directory
        logger.info("Moving mutations and expressions file")
        try:
            shutil.move(mutations_file, self.paths.mutations_file)
            shutil.move(expressions_file, self.paths.expressions_file)
        except PermissionError as e:
            raise FrontendError from e
        except FileNotFoundError as e:
            raise FrontendError from e
        except Exception as e:
            raise BackendError from e
        self.disease = disease
        self.location = location

    def read_mutations_from_file(self,
                                 mut_class=Mutations):
        logger.info("Reading mutations from file")
        get_mutations = mut_class(self.paths.mutations_file)
        mutations = get_mutations.execute()
        return mutations

    def read_expressions_from_file(self,
                                   exp_class=Expressions):
        logger.info("Reading expressions from file")
        get_expressions = exp_class(self.paths.expressions_file)
        expressions = get_expressions.execute()
        return expressions

    def identify_lofs_and_gofs(self,
                               mutations,
                               im_class=IdentifyMutations):
        logger.info("Identifying lof and gof mutations using the mutations file")
        identify_mutations = im_class(mutations, self.paths)
        lofs: set[str]
        gofs: set[str]
        gof_rep: set[str]
        lofs, gofs, gof_rep = identify_mutations.execute()
        return lofs, gofs, gof_rep

    @staticmethod
    def calculate_network_features(expressions,
                                   lofs: set[str],
                                   gofs: set[str],
                                   cnf_class=CalculateNetworkFeatures):
        logger.info("Calculating network features of the gene-gene interaction network")
        calculate_network_features = cnf_class(expressions, lofs, gofs)
        network_features: pd.DataFrame = calculate_network_features.execute()
        return network_features

    @staticmethod
    def predict_essential_genes(network_features,
                                p_class=Predict):
        logger.info("Predicting which genes are essential")
        predictor = p_class(network_features)
        predicted_essential_genes: pd.DataFrame = predictor.execute()
        return predicted_essential_genes

    @staticmethod
    def find_associated_drugs(predicted_essential_genes,
                              gtd_class=GenesToDrugs):
        logger.info("Linking essential genes to cancer drugs")
        genes_to_drugs = gtd_class(predicted_essential_genes)
        essential_drugs_df: pd.DataFrame = genes_to_drugs.execute()
        return essential_drugs_df

    def get_clinical_trials(self,
                            essential_drugs_df):
        logger.info("Searching for relevant clinical trials")
        get_clinical_trials = GetClinicalTrials(essential_drugs_df, self.disease, self.location)
        clinical_trials: dict = get_clinical_trials.execute()
        return clinical_trials

    def get_summary(self, gof_rep: list, lofs: list, essential_drugs_df, clinical_trials):
        logger.info("Summarising information")
        summary = Summary(gof_rep, lofs, essential_drugs_df, clinical_trials)
        summary.write_results_to_file(self.paths.results)
        return summary
