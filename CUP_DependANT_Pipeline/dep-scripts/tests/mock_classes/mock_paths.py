from pathlib import Path

from app.tools.file_paths import Paths, PathException, StaticPaths
from app.tools.config_parser import config


class TestPaths(Paths):
    this_repo = config.this_repo

    data_path = this_repo / 'data'
    if not data_path.exists():
        data_path.mkdir()
    data_dir = str(data_path)

    vep_dir = ""

    @staticmethod
    def find_vep():
        raise NotImplementedError("Should not be using ensembl_vep while testing")

    @staticmethod
    def data(example: str) -> str:
        path = Path(Paths.data_dir) / example
        return str(path)

    @staticmethod
    def reference(example: str) -> str:
        path = Paths.this_repo / 'reference' / example
        if not path.exists():
            raise PathException("Could not find " + str(path))
        return str(path)

    @staticmethod
    def get_latest_classifier() -> Path:
        classifiers_dir = TestPaths.this_repo / "classifiers"

        if config.classifier:
            config_classifier = classifiers_dir / config.classifier
            if config_classifier.exists():
                return config_classifier

        if not classifiers_dir.exists():
            raise PathException(f"Classifiers directory not found: {classifiers_dir}")

        classifier_files = [f for f in classifiers_dir.iterdir() if f.is_file()]

        if not classifier_files:
            raise PathException("No classifier files found in the classifiers directory.")

        latest_file = max(classifier_files, key=lambda f: f.stat().st_mtime)
        return latest_file

    def __init__(self):
        super().__init__("fake_task_id",
                         "mutations_file",
                         "expressions_file")

        def make_patient_path(filename: str) -> str:
            patient_path = str(TestPaths.data_path / filename)
            return patient_path

        self.mutations_file = make_patient_path("mutations.tsv")
        self.expressions_file = make_patient_path("expressions.tsv")
        self.bad_mutations = make_patient_path("bad_mutations.tsv")
        self.results = make_patient_path("results.json")
        self.vep_results = make_patient_path("vep_results.txt")


class TestStaticPaths(StaticPaths):
    oncogenes = TestPaths.reference('oncogenes.txt')
    kinases = TestPaths.reference('kinases.txt')
    tumour_suppressors = TestPaths.reference('tumour_suppressors.txt')
    reactome = TestPaths.reference('reactome.csv')
    chembl = TestPaths.reference('chembl_gene_to_drugs.csv')
    alpha_missense = TestPaths.reference('am_canon.tsv')
    alpha_missense_isoforms = TestPaths.reference('am_isoforms.tsv')
    hugo_conversion = Paths.reference('hugo.json')
