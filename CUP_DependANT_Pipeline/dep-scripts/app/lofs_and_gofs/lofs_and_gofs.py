import pandas as pd

from app.data_pipeline_stage import DataPipelineStage
from app.lofs_and_gofs.ensembl_vep import EnsemblVep
from app.lofs_and_gofs.mutations import Variant, Genes
from app.tools.columns import VepResultDf
from app.tools.file_paths import Paths, StaticPaths
from app.tools.logger import create_logger

logger = create_logger()


class IdentifyMutations(DataPipelineStage):

    def __init__(self, mutations: pd.DataFrame, paths: Paths, ev_class=EnsemblVep, static_paths=StaticPaths):
        logger.debug(f"Initialising IdentifyMutations with patient dir {paths.data_dir}, mutations ~ {mutations.shape}")
        super().__init__(mutations)
        self.paths = paths
        self.genes = Genes(static_paths, mutations)
        self.vep = ev_class(self.paths)

    def execute(self):
        logger.debug("IdentifyMutations.execute()")
        logger.info("Running ensembl-vep")
        vep_results = self.vep.execute()
        logger.info(f"Vep results have {vep_results.shape[0]} rows")

        logger.info("Creating gene object per hugo gene name of ensembl-vep results")
        genes = [self.genes.get_gene(vep_results.loc[vep_results[VepResultDf.SYMBOL] == hugo], hugo)
                 for hugo in vep_results[VepResultDf.SYMBOL].unique() if hugo]
        logger.debug(f"{len(genes)} genes found")

        logger.info("Classifying mutations")
        mutation_type = {Variant.GOF: set(), Variant.LOF: set(), Variant.NONE: set()}
        for gene in genes:
            mutation_type[self.genes.classify(gene)].add(gene)
        logger.debug(f"{len(mutation_type[Variant.GOF])} GOFs, "
                     f"{len(mutation_type[Variant.LOF])} LOFs, "
                     f"{len(mutation_type[Variant.NONE])} other")

        logger.info("Getting lofs")
        lofs = mutation_type[Variant.LOF]
        logger.info("Getting gofs")
        gofs = mutation_type[Variant.GOF]
        logger.info("Getting amino acid substitution from alpha missense where available")
        aa_subs = set()
        for gof in mutation_type[Variant.GOF]:
            aa_subs.update(gof.gofs)

        return {mutation.hugo for mutation in lofs}, {mutation.hugo for mutation in gofs}, aa_subs
