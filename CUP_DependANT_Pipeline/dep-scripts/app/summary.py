import json
from datetime import datetime
import pandas as pd
from app.tools.logger import create_logger
from app.tools.columns import GenesToDrugs

logger = create_logger()


class Summary:

    def __init__(self, gofs: list, lofs: list, drugs: pd.DataFrame, clinical_trials: dict):
        self.date = datetime.now().strftime("%d/%m/%Y")
        self.gofs = gofs
        self.lofs = lofs
        self.drugs = drugs
        self.clinical_trials = clinical_trials

    def get_drug_metadata(self):
        gene_data = {}
        ordered_genes = self.order_genes()
        drug_genes = self.drugs[GenesToDrugs.GENE].unique().tolist()

        groups = self.drugs.groupby(GenesToDrugs.GENE)
        for gene in ordered_genes:
            if gene in drug_genes:
                group = groups.get_group(gene)
                confidence_score = group[GenesToDrugs.CONFIDENCE].iloc[0]
                drugs = group[GenesToDrugs.DRUG].drop_duplicates().tolist()
                gene_data[gene] = {
                    "confidence_score": confidence_score,
                    "drugs": drugs
                }
        return gene_data

    def get_clinical_trials_metadata(self):
        ordered_genes = self.order_genes()
        gene_data = {}
        for gene in ordered_genes:
            if gene in self.clinical_trials.keys() and self.clinical_trials[gene]:
                gene_data[gene] = {}
                for drug in self.clinical_trials[gene].keys():
                    num_trials = len(self.clinical_trials[gene][drug])
                    gene_data[gene][drug] = num_trials
        return gene_data

    def get_results(self):
        summary_dic = {'date': self.date,
                       'gofs': self.gofs,
                       'lofs': self.lofs,
                       'drug_metadata': self.get_drug_metadata(),
                       'clinical_trial_metadata': self.get_clinical_trials_metadata()}
        return summary_dic

    def write_results_to_file(self, filepath):
        summary = self.get_results()
        with open(filepath, 'w') as f:
            json.dump(summary, f, indent=4)
        return summary

    def order_genes(self):
        """ Orders the gene:confidence_score dictionary by confidence_score
            so that higher confidence scores come first
            """
        genes = self.drugs[[GenesToDrugs.GENE, GenesToDrugs.CONFIDENCE]].drop_duplicates()
        essential_genes = genes.set_index('Gene')['confidence_score'].to_dict()
        return dict(list(reversed(sorted(essential_genes.items(), key=lambda x: x[1]))))
