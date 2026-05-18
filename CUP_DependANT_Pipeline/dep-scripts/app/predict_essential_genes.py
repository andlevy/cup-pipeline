import joblib
import pandas as pd

from app.data_pipeline_stage import DataPipelineStage
from app.tools.file_paths import Paths
from app.tools.logger import create_logger

logger = create_logger()


class Predict(DataPipelineStage):

    CUTOFF = 0.5

    def __init__(self, features: pd.DataFrame):
        logger.debug(f"Predict.init(features ~ {features.shape})")
        super().__init__(features)
        self.classifier = joblib.load(Paths.get_latest_classifier())

    def execute(self):

        y = self.classifier.predict_proba(self.data.iloc[:, :12].values)[:, -1]
        predictions = pd.Series(y, index=self.data.index, name="confidence_score")
        predictions = predictions.loc[predictions > Predict.CUTOFF]
        return predictions.to_frame()
