import pandas as pd

from app.data_pipeline_stage import DataPipelineStage
from app.tools.exceptions import IncorrectInputError, BackendError
from app.tools.logger import create_logger

logger = create_logger()


class Expressions(DataPipelineStage):

    def __init__(self, expression_filename: str):
        logger.debug(f"Expressions.__init__({expression_filename})")
        super().__init__(expression_filename)
        self.validate()

    def execute(self):
        logger.debug("Expressions.execute()")
        try:
            expressions = pd.read_csv(self.data,
                                      sep='\t',
                                      header=None,
                                      index_col=0).squeeze("columns")
            return expressions
        except Exception as e:
            logger.exception(f"Change of state of expressions file {self.data} "
                             f"between validation and execution")
            raise BackendError(f"Change of state of expressions file {self.data} "
                               f"between validation and execution") from e

    def validate(self):
        try:
            expressions = pd.read_csv(self.data,
                                      sep='\t',
                                      header=None,
                                      index_col=0).squeeze("columns")
        except Exception as e:
            logger.exception("Failure to read expressions file using pandas")
            raise IncorrectInputError("Failure to read expressions file. "
                                      "Please check that it is correctly spaced for a tsv file") from e
        return expressions
