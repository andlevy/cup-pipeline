import subprocess

import pandas as pd

from app.data_pipeline_stage import DataPipelineStage
from app.tools.columns import VepResultDf, VCF
from app.tools.config_parser import config
from app.tools.exceptions import BackendError
from app.tools.file_paths import Paths
from app.tools.file_reading import read_vcf_as_tsv
from app.tools.logger import create_logger

logger = create_logger()


class EnsemblVep(DataPipelineStage):

    def __init__(self, paths: Paths):
        super().__init__(paths)
        self.vep = Paths.find_vep()

    def execute(self):
        logger.debug("EnsemblVep.execute()")
        result_file = self.run_ensembl_vep()
        return self.process_vep_file(result_file)

    def run_ensembl_vep(self):
        logger.info("Running ensembl-vep from command line with "
                    f"{self.vep} -i {self.data.mutations_file} --force --output_file {self.data.vep_results}"
                    f"--cache --dir {self.data.vep_dir}")
        try:
            s = subprocess.run([self.vep,
                                "-i",
                                self.data.mutations_file,
                                "--force",
                                "--output_file",
                                self.data.vep_results,
                                "--pick",
                                "--pick_order",
                                config.pick_order,
                                "--cache",
                                "--dir",
                                self.data.vep_dir,
                                "--fasta",
                                config.fasta,
                                "--hgvs",
                                "--symbol"],
                               check=False,
                               timeout=600,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        except subprocess.TimeoutExpired as e:
            raise BackendError(f"Ensembl-vep took more than ten minutes to run") from e
        if s.returncode != 0:
            raise BackendError(f"Failed to run ensembl-vep with return code {s.returncode}")
        else:
            logger.info("Ensembl-vep was run with a return code of 0")
            return self.data.vep_results

    def process_vep_file(self, result_file):
        logger.info(f"Reading ensembl-vep results from {result_file}")
        try:
            results = read_vcf_as_tsv(result_file, converters={VepResultDf.CONSEQUENCE: EnsemblVep.parse_set,
                                                               VepResultDf.LOCATION: EnsemblVep.split_location})
        except Exception as e:
            raise BackendError(f"Failed to read VEP results from {result_file}") from e
        results = self.parse_info_column(results)
        results = self.parse_hgvs(results)

        logger.info("Splitting location column into chrom, pos columns")
        results[[VCF.CHROM, VCF.POS]] = pd.DataFrame(results.Location.tolist(), index=results.index)
        pre_filter_size = results.shape[0]
        results = results.drop(results[(results[VCF.CHROM] == "-") | (results[VCF.POS] == "-")].index)
        post_filter_size = results.shape[1]
        if pre_filter_size > 0:
            if post_filter_size == 0:
                logger.error("Cannot parse Ensembl-VEP results. Failed to acquire CHROM and POS columns")
                raise BackendError("Cannot parse Ensembl-VEP results. Failed to acquire CHROM and POS columns")
            if pre_filter_size - post_filter_size > 0:
                logger.warning("Some results are being missed due to incompatible format")
        return results

    def parse_info_column(self, df, info_column_name=VepResultDf.EXTRA):
        logger.info("Breaking down the VEP INFO column into its constituent parts")

        if info_column_name not in df.columns:
            logger.error(f"Vep results file has no column {info_column_name}")
            raise BackendError(f"Vep results file has no column {info_column_name}")

        def parse_info_string(info_string):
            if pd.isna(info_string):
                return {}
            items = info_string.split(';')
            cells = {}
            for item in items:
                if '=' in item:
                    k, v = item.split('=', 1)
                    cells[k] = v
            return cells

        try:
            # Apply the parsing function to create a Series of dictionaries
            info_dicts = df[info_column_name].apply(parse_info_string)
        except Exception as e:
            logger.exception(f"Failed to parse {info_column_name} in Vep results file")
            raise BackendError(f"Failed to parse {info_column_name} in Vep results file") from e

        # Extract all unique keys from the dictionaries
        all_keys = set()
        for info_dict in info_dicts:
            all_keys.update(info_dict.keys())

        if VepResultDf.SYMBOL not in all_keys:
            raise BackendError(f"No HUGO gene ids found in vep results. "
                               f"Please inspect {self.data.vep_results} for clarity")

        if VepResultDf.HGVSC not in all_keys:
            raise BackendError(f"No HGVSc values found in vep results. "
                               f"Please inspect {self.data.vep_results} for clarity")

        # Create new columns and populate them with values from the dictionaries
        for key in all_keys:
            df[key] = info_dicts.apply(lambda d: d.get(key))

        return df

    @staticmethod
    def parse_hgvs(df, hgvs_column=VepResultDf.HGVSC):
        logger.info("Extracting versioned transcript IDs from the HGVSc key in the info column")
        if hgvs_column not in df.columns:
            raise BackendError("No HGVSc values found in vep results. Please inspect vep cli command")

        def parse_row(hgvs_string):
            if pd.isna(hgvs_string):
                return hgvs_string
            return hgvs_string.split(":")[0]

        df[VepResultDf.HGVS] = df.apply(lambda x: parse_row(x[hgvs_column]), axis=1)
        return df

    @staticmethod
    def parse_set(string):
        if isinstance(string, str):
            return set(x.strip() for x in string.split(","))
        return set()

    @staticmethod
    def split_location(loc_str):
        try:
            chromosome, loc = loc_str.split(':')
            return chromosome, int(loc.split("-")[0])
        except ValueError:
            # normally we would throw an exception but we can't because this is a converter method
            # this means we need to return something that is clearly incorrect so it can be picked up later
            logger.warning(f"Could not split location {loc_str}")
            return "-", "-"
