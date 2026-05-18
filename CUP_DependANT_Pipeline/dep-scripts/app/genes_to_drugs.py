from app.data_pipeline_stage import DataPipelineStage
import pandas as pd

from app.tools.file_paths import StaticPaths
from app.tools.logger import create_logger
from app.tools.columns import Chembl

logger = create_logger()


class GenesToDrugs(DataPipelineStage):
    gene_drugs = pd.read_csv(StaticPaths.chembl, index_col=0)

    def __init__(self, essential_genes: pd.DataFrame):
        super().__init__(essential_genes)

    def execute(self):
        drugs = pd.merge(GenesToDrugs.gene_drugs,
                         self.data,
                         right_index=True,
                         left_on=Chembl.GENE)
        drugs.drop_duplicates(subset=[Chembl.DRUG, Chembl.GENE, Chembl.DOC], inplace=True)
        drugs.reset_index(inplace=True)
        return drugs

    @staticmethod
    def get_drug_json(drugs: pd.DataFrame):
        drug_groups = pd.Series(drugs.groupby(by=Chembl.DRUG).groups)

        def groupby_gene(drug_row):
            return pd.Series(drug_row.groupby(by=Chembl.GENE).groups).map(lambda x: drugs.loc[x])

        def supporting_info(row):
            if row[Chembl.DOC] is None:
                return f"IC50={row.IC50:.3g}nM no document reference provided"
            else:
                return f"IC50={row.IC50:.3g}nM document reference {row[Chembl.DOC]}"

        dg = drug_groups.map(lambda x: drugs.loc[x]).map(groupby_gene)
        drug_json = {}
        for drug in dg.index:
            drug_json[drug] = {}
            drug_slice = dg.loc[drug]
            drug_json[drug][Chembl.MAX_PHASE] = f"{drug_slice.iloc[0].max_phase.max():.3g}"
            genes = {}
            for gene in drug_slice.index:
                gene_df = drug_slice.loc[gene]
                genes[gene] = {'confidence_score': f"{gene_df.iloc[0]['confidence_score']:.3g}",
                               'supporting_info': list(gene_df.apply(supporting_info, axis=1).values),
                               'ave_IC50nM': gene_df['IC50'].mean(),
                               'example_doc': f"https://doi.org/{gene_df.iloc[0]['supporting_doc_reference']}"
                               }
            drug_json[drug]['Genes'] = genes
        return drug_json
