import pandas as pd
from app.tools.logger import create_logger
from app.tools.exceptions import BackendError
import json


logger = create_logger()


def get_file(path, description, acquisition):
    logger.info(f"Reading {description} from {path}")
    try:
        return acquisition(path)
    except Exception as e:
        logger.exception(f"Failed to retrieve {description} from file")
        raise BackendError(f"Failed to retrieve {description} from file") from e


def tsv_acquisition(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def gene_list_acquisition(path):
    with open(path) as f:
        return set(f.read().split('\n')[1:])


def json_acquisition(path):
    with open(path) as f:
        return json.load(f)
