from pathlib import Path
from app.tools.config_parser import config
from app.tools.logger import create_logger

logger = create_logger()


class PathException(Exception):

    def __init__(self, message):
        self.message = message


class Paths:
    this_repo = config.this_repo
    logger.info(f"This repo has been identified as {this_repo}")

    data_path = this_repo / 'data'
    if not data_path.exists():
        data_path.mkdir()
    data_dir = str(data_path)
    logger.info(f"Data directory location has been identified as {data_dir}")

    vep_dir = config.cache_dir
    if not Path(vep_dir).exists():
        raise PathException(f"Could not find {vep_dir}")
    logger.info(f"Directory location of ensembl vep cache, plugins etc has been identified as {vep_dir}")

    @staticmethod
    def find_vep():
        path = Paths.this_repo.parent / 'ensembl-vep/vep'
        if not path.exists():
            raise PathException("Could not find " + str(path))
        return str(path)

    @staticmethod
    def data(example: str) -> str:
        path = Path(Paths.data_dir) / example
        logger.debug("data({example}) is {path}".format(example=example, path=str(path)))
        return str(path)

    @staticmethod
    def reference(example: str) -> str:
        path = Paths.this_repo / 'reference' / example
        if not path.exists():
            raise PathException("Could not find " + str(path))
        logger.debug("reference({example}) is {path}".format(example=example, path=str(path)))
        return str(path)

    @staticmethod
    def get_latest_classifier() -> Path:
        logger.info("Finding latest classifier")
        classifiers_dir = Paths.this_repo / "classifiers"

        if config.classifier:
            config_classifier = classifiers_dir / config.classifier
            if config_classifier.exists():
                logger.info(f"Loaded classifier as specified in config: {str(config_classifier)}")
                return config_classifier

        if not classifiers_dir.exists():
            raise PathException(f"Classifiers directory not found: {classifiers_dir}")

        classifier_files = [f for f in classifiers_dir.iterdir() if f.is_file()]
        logger.info(f"{len(classifier_files)} classifiers found in {classifiers_dir}")

        if not classifier_files:
            raise PathException("No classifier files found in the classifiers directory.")

        latest_file = max(classifier_files, key=lambda f: f.stat().st_mtime)
        logger.info(f"Latest classifier file is {latest_file}")
        return latest_file

    def __init__(self,
                 task_id: str,
                 mutations_file,
                 expressions_file):
        root = Path(Paths.data(str(task_id)))
        if not root.is_dir():
            root.mkdir(parents=True, exist_ok=True)

        def make_patient_path(filename: str) -> str:
            patient_path = str(root / filename)
            logger.debug(f"make_patient_path({filename}) is {patient_path}")
            return patient_path

        self.mutations_file = make_patient_path(mutations_file)
        self.expressions_file = make_patient_path(expressions_file)
        self.results = make_patient_path("results.json")
        self.clinical_trials = make_patient_path("clinical_trials.json")
        self.vep_results = make_patient_path("vep_results.txt")


class StaticPaths:
    oncogenes = Paths.reference('oncogenes.txt')
    kinases = Paths.reference('kinases.txt')
    tumour_suppressors = Paths.reference('tumour_suppressors.txt')
    reactome = Paths.reference('reactome.csv')
    chembl = Paths.reference('chembl_gene_to_drugs.csv')
    hugo_conversion = Paths.reference('hugo.json')

    logger.info("All reference paths exist")

    Paths.find_vep()
    logger.info("Ensembl-vep has been located")

    Paths.get_latest_classifier()
    logger.info("There is at least one classifier in classifiers directory")
