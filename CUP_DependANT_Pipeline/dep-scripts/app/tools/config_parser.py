import json
from pathlib import Path
from app.tools.exceptions import BackendError
import sys


class Config:

    def __init__(self,
                 logging_level,
                 classifier,
                 alpha_missense_file,
                 alpha_missense_isoforms,
                 pick_order,
                 canonical_info,
                 cache_dir,
                 fasta,
                 this_repo):
        self.logging_level = logging_level
        self.classifier = classifier
        self.alpha_missense_file = alpha_missense_file
        self.alpha_missense_isoforms = alpha_missense_isoforms
        self.pick_order = pick_order
        self.canonical_info = canonical_info
        self.cache_dir = cache_dir
        self.fasta = fasta
        self.this_repo = this_repo
        print(f"Config read: logging level - {self.logging_level}, "
              f"classifier - {self.classifier}, "
              f"cache_dir - {self.cache_dir}, "
              f"pick_order - {self.pick_order}",
              f"canonical_info - {self.canonical_info}"
              f"fasta - {self.fasta}, "
              f"alpha_missense_file - {self.alpha_missense_file}, "
              f"this_repo - {self.this_repo}")
        self.validate()

    @staticmethod
    def from_json(json_object: dict, repo_address: Path):
        logging_level = json_object["logging_level"]
        classifier = json_object["classifier_name"]
        alpha_missense_file = json_object["alpha_missense_file"]
        alpha_missense_isoforms = json_object["alpha_missense_isoforms"]
        pick_order = json_object["pick_order"]
        canonical_info = json_object["canonical_info"]
        cache_dir = json_object["cache_dir"]
        if cache_dir[-1] == "/":
            cache_dir = cache_dir[:-1]
        if "fasta" in json_object:
            fasta = json_object["fasta"]
        else:
            fasta = f"{cache_dir}/homo_sapiens/110_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
        this_repo = repo_address
        return Config(logging_level,
                      classifier,
                      alpha_missense_file,
                      alpha_missense_isoforms,
                      pick_order,
                      canonical_info,
                      cache_dir,
                      fasta,
                      this_repo)

    def validate(self):
        if self.logging_level not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "STDOUT"]:
            raise BackendError(f"logging level {self.logging_level} isn't valid")
        if not Path(self.cache_dir).exists() or not Path(self.cache_dir).is_dir():
            raise BackendError(f"cache directory {self.cache_dir} is not a valid directory")
        if not Path(self.alpha_missense_file).exists() or not Path(self.alpha_missense_file).is_file():
            raise BackendError(f"alpha missense file {self.alpha_missense_file} is not a valid file")
        if not Path(self.alpha_missense_isoforms).exists() or not Path(self.alpha_missense_isoforms).is_file():
            raise BackendError(f"alpha missense isoforms file {self.alpha_missense_isoforms} is not a valid file")


repo = Path(__file__).resolve().parent.parent.parent
config_file = str(repo / "config.json")

with open(config_file, "r") as f:
    try:
        json_config = json.load(f)
    except Exception as e:
        print(f"Could not read config from {config_file}")
        raise e

if 'unittest' in sys.modules.keys():
    config = Config("STDOUT",
                    "",
                    str(repo / "tests/reference/am_canon.tsv"),
                    str(repo / "tests/reference/am_isoforms.tsv"),
                    "mane_select,mane_plus_clinical,canonical,ccds,biotype,rank",
                    False,
                    repo,
                    "",
                    repo / "tests")
else:
    config = Config.from_json(json_config, repo)
