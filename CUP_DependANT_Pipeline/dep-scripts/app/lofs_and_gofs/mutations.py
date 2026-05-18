from enum import Enum

from app.lofs_and_gofs.alpha_missense import AlphaMissense
from app.tools.columns import VepResultDf, VCF
from app.tools.dataframes import get_file, json_acquisition
from app.tools.exceptions import BackendError
from app.tools.logger import create_logger
import pandas as pd

logger = create_logger()


class Variant(Enum):
    LOF = "lof"
    GOF = "gof"
    NONE = ""


class Gene:

    def __init__(self,
                 rows: pd.DataFrame,
                 hugo):
        self.rows = rows
        self.hugo = hugo
        try:
            self.validate()
        except ValueError as e:
            logger.error("Could not initialise Gene", e)
            raise e
        self.gofs = set()
        logger.debug(f"Gene created with hugo: {self.hugo},"
                     f"rows: {self.rows.shape[0]}")

    def validate(self):
        logger.debug("Gene.validate()")
        if not self.hugo:
            raise ValueError("No HUGO symbol provided")
        # 1. All rows should be on the same chromosome
        # 2. shouldn't have different alleles with the same [chrom, position]
        if self.rows[VCF.CHROM].nunique() != 1:
            raise ValueError(f"Mutations associated with {self.hugo} lie on more than one chromosome")
        if self.rows[VCF.POS].nunique() != self.rows.shape[0]:
            raise ValueError(f"Duplicate mutations found in the VCF file")

    @staticmethod
    def concatenate_sets(series):
        all_values = set()
        for a_set in series:
            if isinstance(a_set, set):
                all_values.update(a_set)
            else:
                raise ValueError(f"Unable to parse vep results: consequence {a_set} is not a set")
        return all_values


class Mutation:

    def __init__(self, chromosome, position, consequences, hugo, versioned_transcript):
        self.chromosome = chromosome
        self.position = position
        self.consequences = Mutation.ignore_dash(consequences)
        self.hugo = hugo
        self.versioned_transcript = versioned_transcript
        self.extra = ""
        logger.debug(f"Mutation created with chromosome: {self.chromosome},"
                     f"position: {self.position},"
                     f"consequences: {self.consequences},"
                     f"hugo: {self.hugo}")

    def set_extra(self, extra):
        self.extra = extra

    @staticmethod
    def ignore_dash(attribute):
        if attribute == "-":
            return ""
        return attribute

    def get_aa_sub(self):
        if self.hugo:
            if self.extra:
                return f"{self.hugo}-{self.extra}"
            return self.hugo
        raise ValueError(f"No Hugo gene name for this mutation {self.chromosome}:{self.position}")

    def is_pathogenic(self, mutations: pd.DataFrame):
        logger.debug(f"Mutation.is_pathogenic() for {self.chromosome}:{self.position}")
        if VepResultDf.missense_consequence not in self.consequences:
            raise ValueError("Cannot call Alpha Missense with a mutation that is not a missense variant")
        row = mutations.loc[(mutations[VCF.CHROM] == int(self.chromosome)) &
                            (mutations[VCF.POS] == int(self.position))].iloc[0]

        ref, alt = row[VCF.REF], row[VCF.ALT]
        alpha_missense = AlphaMissense(self.chromosome, self.position, ref, alt, self.versioned_transcript)
        status, extra = alpha_missense.double_query()
        if not status:
            return False
        self.set_extra(extra)
        return True

    @staticmethod
    def create_mutation(row, hugo):
        logger.debug(f"Mutation.create_mutation({row})")
        return Mutation(row[VCF.CHROM],
                        row[VCF.POS],
                        row[VepResultDf.CONSEQUENCE],
                        hugo,
                        row[VepResultDf.HGVS])

    def __eq__(self, other):
        if isinstance(other, Mutation):
            return (self.chromosome == other.chromosome and
                    self.position == other.position and
                    self.hugo == other.hugo)
        return False

    def __hash__(self):
        return hash((self.chromosome, self.position, self.hugo))


class Genes:

    def __init__(self, static_paths, original_vcf: pd.DataFrame):
        try:
            self.oncogenes = set(get_file(static_paths.oncogenes, "list of oncogenes", json_acquisition))
            self.kinases = set(get_file(static_paths.kinases, "list of kinases", json_acquisition))
            self.tumour_suppressors = set(get_file(static_paths.tumour_suppressors,
                                                   "list of tumour_suppressors",
                                                   json_acquisition))
            self.hugo_conversion = get_file(static_paths.hugo_conversion,
                                            "Hugo conversion dictionary",
                                            json_acquisition)
        except BackendError as e:
            logger.error("Couldn't initialise Genes due to failure to load reference files", e)
            raise e
        except Exception as e:
            logger.error("Couldn't initialise Genes due to unknown error", e)
            raise e
        self.original_vcf = original_vcf

    def get_gene(self, rows: pd.DataFrame, gene_name):
        return Gene(rows, self.get_uptodate_hugo(gene_name))

    def get_uptodate_hugo(self, gene_name):
        if gene_name in self.hugo_conversion:
            return self.hugo_conversion[gene_name]
        return None

    @staticmethod
    def is_bla(hugo, gene_list, description):
        if hugo:
            if hugo in gene_list:
                logger.debug(f"{hugo} found in {description} list")
                return True
        return False

    def is_oncogene(self, hugo):
        return Genes.is_bla(hugo, self.oncogenes, "oncogenes")

    def is_kinase(self, hugo):
        return Genes.is_bla(hugo, self.kinases, "kinases")

    def is_tumour_suppressor(self, hugo):
        return Genes.is_bla(hugo, self.tumour_suppressors, "tumour_suppressors")

    # How a pathogenic missense mutation from this gene would be classified
    def classify_pathogenic_missense(self, hugo):
        if self.is_tumour_suppressor(hugo):
            return Variant.LOF
        if self.is_oncogene(hugo) or self.is_kinase(hugo):
            return Variant.GOF
        return Variant.LOF

    def classify(self, gene: Gene):
        logger.debug(f"Genes.classify of {gene.hugo}")
        if not gene.hugo:
            logger.warning(
                f"Ignoring {gene.rows.shape[0]} rows because they do not have a depmap compatible hugo name.")
            return Variant.NONE
        consequences = Gene.concatenate_sets(gene.rows[VepResultDf.CONSEQUENCE])
        if len(VepResultDf.lof_consequences.intersection(consequences)) != 0:
            logger.info(f"{gene.hugo} has a LOF mutation")
            return Variant.LOF

        logger.info("Creating mutation object per row of ensembl-vep results")
        mutations = [Mutation.create_mutation(row, gene.hugo) for idx, row in gene.rows.iterrows()]
        for mutation in mutations:
            if VepResultDf.missense_consequence in mutation.consequences:
                if mutation.is_pathogenic(self.original_vcf):
                    variant = self.classify_pathogenic_missense(gene.hugo)
                    if variant == Variant.LOF:
                        logger.info(f"{gene.hugo} has a missense LOF mutation")
                        return Variant.LOF
                    gene.gofs.add(mutation.get_aa_sub())
        if len(gene.gofs) > 0:
            logger.info(f"{gene.hugo} has a GOF mutation")
            return Variant.GOF
        return Variant.NONE
