#!/usr/bin/env bash
#
# download_open_targets.sh
#
# Downloads selected Open Targets Platform datasets (Parquet format) from the
# EBI FTP mirror, into a directory you control.
#
# Open Targets release the full Platform as a set of separate dataset
# directories, e.g.:
#   https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/<RELEASE>/output/<DATASET>/
#
# Usage:
#   ./download_open_targets.sh -o /path/to/data
#   ./download_open_targets.sh -o /path/to/data -r 25.06 -d target -d disease
#   ./download_open_targets.sh -r 25.06 --list      # see what's available first
#
# Requires: wget, curl
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults — override with flags below
# ---------------------------------------------------------------------------
RELEASE="26.06"
OUTPUT_DIR="./opentargets_data"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
LIST_ONLY=0

# These 15 are the set chosen after reviewing the full Croissant metadata for
# release 26.06 (see README.md in this directory for what each one provides).
# Override with -d to download a different set instead.
DEFAULT_DATASETS=(
  target
  disease
  target_essentiality
  target_prioritisation
  drug_mechanism_of_action
  drug_molecule
  clinical_target
  clinical_indication
  clinical_report
  evidence_cancer_gene_census
  evidence_intogen
  evidence_cancer_biomarkers
  evidence_crispr_screen
  evidence_eva_somatic
  association_by_datasource_direct
)
DATASETS=()

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [-o OUTPUT_DIR] [-r RELEASE] [-d DATASET ...] [--list]

  -o, --output DIR     Where to save the data (default: ${OUTPUT_DIR})
  -r, --release VER    Open Targets release version, e.g. 26.06
                        (default: ${RELEASE}; check ${BASE_URL}/
                        for newer releases)
  -d, --dataset NAME    Dataset to download. Repeatable. If omitted, the
                        default set is used (the 15 datasets decided on for
                        this project):
                        ${DEFAULT_DATASETS[*]}
  -l, --list            List datasets available for the chosen release and exit
  -h, --help            Show this help

Examples:
  $0 --list -r 26.06
  $0 -o ~/data/opentargets                  # downloads the default 15
  $0 -o ~/data/opentargets -d target -d disease
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
    -r|--release) RELEASE="$2"; shift 2 ;;
    -d|--dataset) DATASETS+=("$2"); shift 2 ;;
    -l|--list) LIST_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ ${#DATASETS[@]} -eq 0 ]]; then
  DATASETS=("${DEFAULT_DATASETS[@]}")
fi

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for tool in wget curl; do
  if ! command -v "$tool" &>/dev/null; then
    echo "Required tool '$tool' is not installed. Please install it and re-run." >&2
    exit 1
  fi
done

RELEASE_URL="${BASE_URL}/${RELEASE}/output"

echo "Checking release ${RELEASE} exists at ${RELEASE_URL} ..."
if ! curl -sf --max-time 15 -o /dev/null "${RELEASE_URL}/"; then
  echo "Could not reach ${RELEASE_URL}/ — check the release version and your network access." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# --list: show available dataset directories for this release, then exit
# ---------------------------------------------------------------------------
if [[ "$LIST_ONLY" -eq 1 ]]; then
  echo "Available datasets for release ${RELEASE}:"
  curl -s --max-time 30 "${RELEASE_URL}/" \
    | grep -oP 'href="([^"/?]+)/"' \
    | sed -E 's/href="(.+)\/"/\1/' \
    | sort -u
  exit 0
fi

# ---------------------------------------------------------------------------
# Download each requested dataset
# ---------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
echo "Downloading release ${RELEASE} datasets into: ${OUTPUT_DIR}"
echo "Datasets: ${DATASETS[*]}"
echo

FAILED=()
for dataset in "${DATASETS[@]}"; do
  DATASET_URL="${RELEASE_URL}/${dataset}"
  echo "----------------------------------------------------------------"
  echo "Fetching: ${dataset}"
  echo "From:     ${DATASET_URL}/"
  echo "To:       ${OUTPUT_DIR}/${dataset}/"

  if wget --quiet --show-progress \
      --recursive --no-parent --no-host-directories \
      --cut-dirs 8 \
      --directory-prefix="${OUTPUT_DIR}/${dataset}" \
      "${DATASET_URL}/"; then
    echo "Done: ${dataset}"
  else
    echo "FAILED: ${dataset}" >&2
    FAILED+=("$dataset")
  fi
done

echo
echo "================================================================"
if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "All datasets downloaded successfully to: ${OUTPUT_DIR}"
else
  echo "Completed with failures in: ${FAILED[*]}"
  echo "Check the dataset name is correct for release ${RELEASE} (use --list)."
  exit 1
fi
