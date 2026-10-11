#!/usr/bin/env bash
# =============================================================
# 00_download_data.sh
#
# Download the raw 10x Genomics matrices for GSE155698
# (Peng et al., PDAC single-cell RNA-seq) and lay them out in the
# structure expected by scripts/01_load_and_merge_samples.R:
#
#     data/raw/<sample>/filtered_feature_bc_matrix/
#
# Only the ten samples used by this analysis are fetched
# (7 tumour + 3 adjacent-normal). PBMC samples in GSE155698 are
# not used here.
#
# Usage (from the repository root):
#     bash scripts/00_download_data.sh
# =============================================================

set -euo pipefail

RAW_DIR="data/raw"
mkdir -p "$RAW_DIR"

# GSM accession -> sample name used in data/sample_metadata.csv
declare -A SAMPLES=(
  [GSM4710691]=PDAC_TISSUE_3
  [GSM4710693]=PDAC_TISSUE_5
  [GSM4710698]=PDAC_TISSUE_10
  [GSM4710699]=PDAC_TISSUE_11A
  [GSM4710701]=PDAC_TISSUE_12
  [GSM4710704]=PDAC_TISSUE_15
  [GSM4710705]=PDAC_TISSUE_16
  [GSM4710706]=AdjNorm_TISSUE_1
  [GSM4710707]=AdjNorm_TISSUE_2
  [GSM4710708]=AdjNorm_TISSUE_3
)

for gsm in "${!SAMPLES[@]}"; do
  sample="${SAMPLES[$gsm]}"
  prefix="GSM4710nnn/${gsm}"
  ssh_stem="${gsm}_${sample}"
  url="https://ftp.ncbi.nlm.nih.gov/geo/samples/${prefix}/suppl/${ssh_stem}.tar.gz"

  if [[ -d "${RAW_DIR}/${sample}" ]]; then
    echo "[skip] ${sample} already present"
    continue
  fi

  echo "[download] ${sample} (${gsm})"
  curl -4 -sSL -C - -o "${RAW_DIR}/${ssh_stem}.tar.gz" "${url}"

  # Extract; the archive contains a top-level directory named after the
  # sample (e.g. PDAC_TISSUE_3/).
  tar xzf "${RAW_DIR}/${ssh_stem}.tar.gz" -C "${RAW_DIR}"
  rm -f "${RAW_DIR}/${ssh_stem}.tar.gz"
done

echo
echo "Downloaded samples:"
ls -1 "$RAW_DIR"
