import pandas as pd

from app.tools.logger import create_logger

logger = create_logger()


def read_vcf_as_tsv(filename, **kwargs):
    skip_rows = 0
    with open(filename, "r+") as f:
        while next(f)[0:2] == "##":
            skip_rows += 1
    df = pd.read_csv(filename, sep="\t", skiprows=skip_rows, **kwargs)
    return df
