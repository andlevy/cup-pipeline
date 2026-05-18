import pandas as pd

from app.clinical_trials.clinical_trials import find_all_clinical_trials
from app.data_pipeline_stage import DataPipelineStage
from app.tools.logger import clinical_trials_logger
from app.tools.exceptions import IncorrectInputError, FrontendError

logger = clinical_trials_logger()


class GetClinicalTrials(DataPipelineStage):

    def __init__(self, drugs: pd.DataFrame, disease: str, location: str):
        logger.debug(f"Initialising GetClinicalTrials with {drugs}, {disease}, {location}")
        super().__init__(drugs)
        self.disease = disease
        self.location = location

    def __str__(self):
        return f"GetClinicalTrials: disease={self.disease}, location={self.location}"

    def __repr__(self):
        return f"GetClinicalTrials({self.disease!r}, {self.location!r})"

    def execute(self):
        logger.debug("GetClinicalTrials.execute()")
        drugs = self.data["Drug"].unique()
        drug_trials_df = []
        trials_by_gene_by_drug = {}
        for drug in drugs:
            logger.debug(f"Finding clinical trials for {drug}")
            genes = self.data.loc[self.data["Drug"] == drug].drop_duplicates()
            gene = genes[genes['confidence_score'] == genes['confidence_score'].max()]['Gene'].iloc[0]
            trials = find_all_clinical_trials(self.disease,
                                              [drug],
                                              gene,
                                              self.location)
            if len(trials) != 0:
                if gene not in trials_by_gene_by_drug.keys():
                    trials_by_gene_by_drug[gene] = {drug: [trial.create_dict() for trial in trials]}
                else:
                    trials_by_gene_by_drug[gene][drug] = [trial.create_dict() for trial in trials]
                drug_trials_df += [[gene,
                                    drug,
                                    trial.title,
                                    trial.url,
                                    trial.summary] for trial in trials]
        return trials_by_gene_by_drug

    def validate(self):
        if not isinstance(self.disease, str):
            raise FrontendError("disease not of type string")
        if "%" in self.disease:
            raise IncorrectInputError("Invalid character % in disease string")
        if not isinstance(self.location, str):
            raise FrontendError("location not of type string")
        if "%" in self.location:
            raise IncorrectInputError("Invalid character % in location string")
