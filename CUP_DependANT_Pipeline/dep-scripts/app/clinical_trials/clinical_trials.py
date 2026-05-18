from app.clinical_trials.eligibility import format_eligibility
import requests
from app.tools.logger import clinical_trials_logger

logger = clinical_trials_logger()


class ClinicalTrials:

    def __init__(self,
                 url: str,
                 drugs: list[str],
                 gene: str,
                 title: str,
                 summary: str,
                 conditions: list[str],
                 countries: list[str],
                 eligibility: list[str]):
        self.url = url
        self.drugs = ' '.join(drugs)
        self.gene = gene,
        self.title = title
        self.summary = summary
        self.conditions = conditions
        self.countries = countries
        self.eligibility = eligibility

    def __str__(self):
        return (f"ClinicalTrials: url={self.url}, drugs={self.drugs}, gene={self.gene}, title={self.title}, "
                f"summary={self.summary}, conditions={self.conditions}, countries={self.countries}, "
                f"eligibility={self.eligibility[0]} ...")

    def __repr__(self):
        return (f"ClinicalTrials: {self.url!r}, {self.drugs!r}, {self.gene!r}, {self.title!r}, "
                f"{self.summary!r}, {self.conditions!r}, {self.countries!r}, "
                f"{self.eligibility[0]!r} ...")

    def create_dict(self):
        logger.debug("clinical_trials.create_dict()")
        self_dict = {"url": self.url,
                     "drug": self.drugs,
                     "gene": self.gene,
                     "title": self.title,
                     "summary": self.summary,
                     "conditions": self.conditions,
                     "countries": self.countries,
                     "eligibility": self.eligibility}
        return self_dict

    @staticmethod
    def create_url_from_nct_id(nct_id: str):
        logger.debug(f"clinical_trials.create_url_from_nct_id({nct_id})")
        return "https://clinicaltrials.gov/study/" + nct_id

    @staticmethod
    def create_with_json(json: dict, drugs: list[str], gene: str):
        logger.debug(f"clinical_trials.create_with_json(<json>, {drugs}, {gene})")
        if not ClinicalTrials.valid_clinical_trial_json(json):
            return None
        protocol_section = json['protocolSection']
        derived_section = json['derivedSection']
        if 'eligibilityCriteria' in protocol_section['eligibilityModule']:
            eligibility = format_eligibility(protocol_section['eligibilityModule']['eligibilityCriteria'])
        else:
            eligibility = []
        nct_id = json['protocolSection']['identificationModule']['nctId']
        if 'briefSummary' in protocol_section['descriptionModule']:
            summary = protocol_section['descriptionModule']['briefSummary']
        else:
            summary = ""
        countries = list(set([i['country'] for i in protocol_section['contactsLocationsModule']['locations']]))
        conditions = list(set([i["asFound"] for i in derived_section['conditionBrowseModule']['browseLeaves'] if
                               i["relevance"] == "HIGH"]))
        if 'briefTitle' in protocol_section['identificationModule']:
            title = protocol_section['identificationModule']['briefTitle']
        elif 'title' in protocol_section['identificationModule']:
            title = protocol_section['identificationModule']['title']
        else:
            title = ""
        return ClinicalTrials(ClinicalTrials.create_url_from_nct_id(nct_id),
                              drugs,
                              gene,
                              title,
                              summary,
                              conditions,
                              countries,
                              eligibility)

    @staticmethod
    def format_eligibility(eligibility_string: str):
        logger.debug(f"clinical_trials.format_eligibility({eligibility_string})")
        return [line for line in eligibility_string.splitlines()]

    @staticmethod
    def valid_clinical_trial_json(json: dict):
        logger.debug(f"clinical_trials.valid_clinical_trial_json(<json>)")
        if 'protocolSection' not in json or 'derivedSection' not in json:
            logger.debug("clinical trials json is invalid - missing protocol_Section or derivedSection")
            return False
        protocol_section = json['protocolSection']
        derived_section = json['derivedSection']
        if 'eligibilityModule' in protocol_section and \
                'descriptionModule' in protocol_section and \
                'identificationModule' in protocol_section and \
                'contactsLocationsModule' in protocol_section and \
                'conditionBrowseModule' in derived_section:
            if 'nctId' in protocol_section['identificationModule'] and \
                    'locations' in protocol_section['contactsLocationsModule'] and \
                    'browseLeaves' in derived_section['conditionBrowseModule']:
                logger.debug("clinical trials json is valid")
                return True
        logger.debug("clinical trials json is invalid - missing submodules")
        return False


def clinical_trials_request(disease: str,
                            drugs: list,
                            location: str,
                            trial_status=["RECRUITING"],
                            page_token=""):
    logger.debug(f"clinical_trials.clinical_trials_request({disease}, "
                 f"{drugs}, {location}, {trial_status}, {page_token})")
    api_string = "https://clinicaltrials.gov/api/v2/studies?"
    if disease is not None and disease != "":
        api_string += "query.cond=%22{disease}%22".format(disease=disease.replace(" ", "+"))
    if drugs is not None and drugs:
        drug_string = ' '.join(drugs)
        api_string += "&query.term={drugs}".format(drugs=drug_string.replace(" ", "+"))
    if location not in [None, "", "All"]:
        api_string += "&query.locn=%22{location}%22".format(location=location.replace(" ", "+"))
    if trial_status is not None and trial_status:
        status_string = "%2C+%7C".join(trial_status)
        api_string += "&filter.overallStatus={trial_status}".format(trial_status=status_string)
    if page_token is not None and page_token:
        api_string += "&pageToken={token}".format(token=page_token)
    logger.info(f"Requesting {api_string}")
    response = requests.get(api_string)
    if response.status_code != 200:
        logger.warning(f"Unexpected status code from {api_string}")
        return None
    if response.headers["content-type"].strip() != 'application/json':
        logger.warning(f"Unexpected content type from {api_string}")
        return None
    return response.json()


def retrieve_clinical_trials_from_response(data, drugs, gene):
    logger.debug(f"clinical_trials.retrieve_clinical_trials_from_response({data}, {drugs}, {gene})")
    clinical_trials = []
    if 'studies' in data and len(data['studies']) > 0:
        for study in data['studies']:
            if ClinicalTrials.valid_clinical_trial_json(study):
                clinical_trials.append(ClinicalTrials.create_with_json(study, drugs, gene))
    else:
        if 'studies' not in data:
            logger.warning("Could not find studies from JSON")
        else:
            logger.debug("Studies was empty")
    return clinical_trials


def find_all_clinical_trials(disease: str,
                             drugs: list,
                             gene: str,
                             location: str,
                             trial_status=["RECRUITING"]):
    logger.debug(f"clinical_trials.find_all_clinical_trials({disease}, {drugs}, {gene}, {location}, {trial_status})")
    data = clinical_trials_request(disease=disease,
                                   drugs=drugs,
                                   location=location,
                                   trial_status=trial_status)
    clinical_trials = retrieve_clinical_trials_from_response(data, drugs, gene)
    i = 0
    while 'nextPageToken' in data and i < 10:
        i += 1
        logger.debug(f"Attempting to retrieve more data using nextPageToken {data['nextPageToken']}")
        data = clinical_trials_request(disease=disease,
                                       drugs=drugs,
                                       location=location,
                                       trial_status=trial_status,
                                       page_token=data['nextPageToken'])
        new_clinical_trials = retrieve_clinical_trials_from_response(data, drugs, gene)
        logger.debug(f"{len(new_clinical_trials)} retrieved")
        clinical_trials += new_clinical_trials
    logger.info(f"{len(clinical_trials)} found for {disease}, {drugs}, {gene}, {location}, {trial_status}")
    return clinical_trials
