# CUP Gene Target Database

Builds a local copy of selected [Open Targets Platform](https://platform.opentargets.org/) datasets. These are used to prioritise candidate drug targets from CUP sequencing results.

## Usage

```bash
conda env create -f ../envs/OTP_GeneSelection.yml   # pandas, pyarrow, pyspark

./download_open_targets.sh --list -r 26.06          # list available datasets
./download_open_targets.sh -o /path/to/opentargets  # download the default 15 datasets
./download_open_targets.sh -o /path/to/opentargets -d target -d disease   # custom set
```

Data is downloaded as Parquet from the EBI FTP mirror:
`https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/<RELEASE>/output/<DATASET>/`

## Default Datasets (release 26.06)

| Group | Dataset | Content |
|-------|---------|---------|
| Core entities | `target` | Gene/target annotation (Ensembl IDs, symbols, tractability, location) |
| | `disease` | Disease/phenotype ontology (EFO) |
| Target biology | `target_essentiality` | Gene essentiality from CRISPR screens (DepMap) |
| | `target_prioritisation` | Target prioritisation factors (safety, tractability, etc.) |
| Drugs | `drug_molecule` | Drug/molecule records (ChEMBL) |
| | `drug_mechanism_of_action` | Drug → target mechanism of action |
| Clinical | `clinical_target` | Clinical-stage evidence by target |
| | `clinical_indication` | Clinical-stage evidence by indication |
| | `clinical_report` | Clinical trial / report records |
| Cancer evidence | `evidence_cancer_gene_census` | COSMIC Cancer Gene Census evidence |
| | `evidence_intogen` | IntOGen cancer driver evidence |
| | `evidence_cancer_biomarkers` | Cancer Genome Interpreter biomarkers |
| | `evidence_crispr_screen` | Cancer CRISPR screen evidence |
| | `evidence_eva_somatic` | ClinVar somatic (EVA) evidence |
| Associations | `association_by_datasource_direct` | Direct target–disease association scores per data source |

The set was chosen after reviewing the release's Croissant metadata. Datasets focused on germline genetics, literature mining and non-cancer evidence were left out.
