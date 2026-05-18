from app.data_pipeline_stage import DataPipelineStage
from math import tanh, log
import pandas as pd
from app.gene_gene_interaction import GGI
from app.tools.file_paths import StaticPaths
from app.tools.logger import create_logger
from app.tools.dataframes import get_file, json_acquisition
from app.tools.exceptions import BackendError, IncorrectInputError

logger = create_logger()


class CalculateNetworkFeatures(DataPipelineStage):

    @staticmethod
    def convert_expression_genes(expressions: pd.Series, hugo_dict: dict[str, str]):
        """
        We need the gene names that the expressions are mapped against to match the gene names of the mutations
        and of our gene-gene interaction data
        :param expressions: Indexed by gene name, value is expression
        :param hugo_dict: Keys are old, new and aliases of hugo gene names, values are current hugo gene names
        :return: expressions indexed by the most current hugo gene name
        """
        try:
            conversion = expressions.index.map(hugo_dict)
            relevant_expressions = expressions.loc[conversion.notna()]
            if relevant_expressions.shape[0] == 0:
                raise IncorrectInputError("No matches found between the gene names in the expressions file and our"
                                          "reference list of gene aliases. Please check your file and consider updating"
                                          "to the most recent HGNC version for your gene names.")
            relevant_expressions.index = conversion.dropna()
            return relevant_expressions
        except IncorrectInputError:
            raise
        except Exception as e:
            raise BackendError("Failed to convert expressions file to most recent version of HGNC") from e

    def __init__(self,
                 expressions: pd.Series,
                 lof: set[str],
                 gof: set[str],
                 static_paths=StaticPaths):
        logger.debug(f"CalculateNetworkFeatures.init(expressions ~ {expressions.shape}, "
                     f"lof ~ {len(lof)}, gof ~ {len(gof)})")
        hugo_conversion = get_file(static_paths.hugo_conversion, "hugo conversion dict", json_acquisition)
        if not hugo_conversion or not isinstance(hugo_conversion, dict):
            raise BackendError(f"Hugo conversion data from {static_paths.hugo_conversion} "
                               f"not in correct format (dict[str, str])")
        super().__init__(CalculateNetworkFeatures.convert_expression_genes(expressions, hugo_conversion))
        self.lof = lof
        self.gof = gof
        if len(self.lof.intersection(self.gof)) != 0:
            raise BackendError("Unexpected behaviour: received genes that are marked lof and gof")

    def get_scores(self) -> pd.DataFrame:
        logger.debug("CalculateNetworkFeatures._get_scores()")
        """
        Bring together information from gof lof and expressions
        The higher the expression the lower the score, zeros are reserved for lofs which will be removed later.
        """

        def score(x):
            if x == 0:
                return 0
            else:
                return 0.5 - 0.5 * (tanh(log(x)))

        t = pd.concat([pd.Series(0.0, index=list(self.lof)),
                       pd.Series(10.0, index=list(self.gof)),
                       pd.Series(1.0, index=list(set(self.data.index) - set(self.gof) - set(self.lof)))], axis=0)

        adjusted = t[self.data.index] * self.data
        return adjusted.map(score)

    def execute(self):
        logger.debug("CalculateNetworkFeatures.execute()")
        try:
            scores = self.get_scores()
            ggi = GGI(scores)
            return ggi.get_features()
        except Exception as e:
            logger.exception("Failed to calculate network features")
            raise BackendError("Failed to calculate network features") from e
