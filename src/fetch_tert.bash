#!/usr/bin/env bash
set -euo pipefail

# Resolve paths relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"

# Output files
out_nt="${INPUT_DIR}/plant_TERT_regions.fasta"
out_aa="${INPUT_DIR}/plant_TERT_regions.faa"

# If both outputs already exist and are non-empty, skip unless FORCE=1
if [[ "${FORCE:-0}" != "1" ]] && [[ -s "$out_nt" && -s "$out_aa" ]]; then
  echo "[fetch_tert] Outputs exist ($out_nt, $out_aa); skipping. Set FORCE=1 to regenerate."
  exit 0
fi

echo "[fetch_tert] Writing nucleotide regions to $out_nt"
: > "$out_nt"

# Define each sequence with: accession, start, end, strand, and name
entries=(
  "NC_003076.8 5538180 5543657 - Arabidopsis_thaliana"
  "NC_050104.1 134209849 134227624 + Zea_mays"
  "NC_089046.1 11424226 11440551 - Oryza_sativa"
  "NW_015810241.1 1616 13759 + Nicotiana_tabacum"
  "NC_037253.2 30018133 30025820 + Physcomitrium_patens"
  "NC_091065.1 117563160 117575366 - Nicotiana_sylvestris"
  "NC_090812.1 821889 841870 - Nicotiana_tomentosiformis"
  "NC_068337.1 10074968 10091244 - Oryza_glaberrima"
  "NW_003314277.1 1528125 1532168 - Selaginella_moellendorffii"
  "NW_019014448.1 285542 300175 + Carica_papaya"
  "NW_021319083.1 5746516 5797768 - Dendrobium_catenatum"
  "NC_037287.2 2358001 2367828 - Populus_trichocarpa"
)

for entry in "${entries[@]}"; do
  read -r acc start end strand species <<< "$entry"

  strand_flag=1
  [[ "$strand" == "-" ]] && strand_flag=2

  /software/team301/edirect/efetch \
    -db nuccore -id "$acc" -format fasta \
    -seq_start "$start" -seq_stop "$end" -strand "$strand_flag" \
    | sed "s/^>/>${species}_/" >> "$out_nt"
done

echo "[fetch_tert] Saved nucleotide regions to $out_nt"

echo "[fetch_tert] Fetching amino acid TERTs to $out_aa"
/software/team301/edirect/esearch \
    -db gene -query 'TERT[Gene Name] AND "Viridiplantae"[Organism]' \
  | /software/team301/edirect/elink -target protein \
  | /software/team301/edirect/efetch -format fasta > "$out_aa"

echo "[fetch_tert] Done."

