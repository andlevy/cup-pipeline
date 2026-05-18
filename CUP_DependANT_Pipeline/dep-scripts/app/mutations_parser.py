from app.data_pipeline_stage import DataPipelineStage
from app.tools.file_reading import read_vcf_as_tsv
from app.tools.columns import VCF
from app.tools.exceptions import IncorrectInputError, BackendError
from app.tools.logger import create_logger

logger = create_logger()


class Mutations(DataPipelineStage):

    def __init__(self, mutation_filename: str):
        super().__init__(mutation_filename)
        self.validate()

    def execute(self):
        logger.debug("Mutations.execute()")
        try:
            mutations = read_vcf_as_tsv(self.data)
            return mutations[[VCF.CHROM, VCF.POS, VCF.REF, VCF.ALT]]
        except Exception as e:
            logger.exception(f"Change of state of mutations file {self.data} "
                             f"between validation and execution")
            # We've already validated so something is really wrong if we get an error now
            raise BackendError(f"Change of state of mutations file {self.data} "
                               f"between validation and execution") from e

    def validate(self):
        try:
            mutations = read_vcf_as_tsv(self.data)
        except Exception as e:
            logger.exception(f"Failure to read mutations file {self.data} using pandas")
            raise IncorrectInputError(f"Failure to read mutations file. "
                                      f"Please check that file is correctly spaced for a tsv file.") from e

        required_columns = [VCF.CHROM,
                            VCF.POS,
                            VCF.REF,
                            VCF.ALT,
                            VCF.ID,
                            VCF.INFO,
                            VCF.QUALITY,
                            VCF.FILTER]
        for i in required_columns:
            if i not in mutations.columns:
                logger.error(f"vcf file not valid: missing column {i}. "
                             "Please provide a tab separated file with"
                             "#CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns")
                raise IncorrectInputError(f"mutations file not valid: missing column {i}. "
                                          "Please provide a tab separated file with"
                                          "#CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns")
        return mutations
