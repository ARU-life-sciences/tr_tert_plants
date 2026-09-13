#!/usr/bin/env bash
set -euo pipefail

# ./src/run_extract_tr_domains.bash
#
# Runs extract_tr_domains.py for every species with an nhmmer TR hit
# (outputs/tr_tbls/<species>_all_tr_seqs_020226.tbl), extracting each TR
# gene domain (USE, TATA-box, G-rich 5' end, Template, Conserved region,
# C-rich 3' end, Terminator - see inputs/TR_domains.all_tr_seqs.tsv, built
# by map_tr_domains.py) from the species' genome assembly, with the
# Template cross-checked against its TIDK-discovered telomeric repeat.
#
# Local (no bsub) - genome decompression dominates runtime, so this runs
# several species in parallel. Safe to re-run: skips species that already
# have an output file.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
# Note: inputs/extra_genomes.txt (Santalum album) is not included here - it
# was only ever run through the TERT tblastn side of the pipeline
# (run_tblastn_tert_local.bash), never through run_tr_hmms.bash, so there is
# no TR nhmmer tbl for it and its genome file lives on a local machine, not
# this cluster.
TBL_DIR="${OUTPUT_DIR}/tr_tbls"
DOMAINS_TSV="${DOMAINS_TSV:-${INPUT_DIR}/TR_domains.all_tr_seqs.tsv}"
TIDK_DIR="${OUTPUT_DIR}/tidk"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_domains}"
JOBS="${JOBS:-12}"

mkdir -p "$OUTDIR"

run_one () {
  species="$1"; fasta="$2"
  out="${OUTDIR}/${species}.tsv"
  if [[ -s "$out" ]]; then
    echo "[run_extract_tr_domains] ${species}: $out exists, skipping."
    return
  fi
  tbl="${TBL_DIR}/${species}_all_tr_seqs_020226.tbl"
  if [[ ! -s "$tbl" ]]; then
    echo "[run_extract_tr_domains] ${species}: no tbl at $tbl, skipping."
    return
  fi
  tidk="${TIDK_DIR}/${species}.tidk.tsv"
  [[ -s "$tidk" ]] || tidk=""
  python3 "${SCRIPT_DIR}/extract_tr_domains.py" "$species" "$fasta" "$tbl" "$DOMAINS_TSV" "$tidk" "$OUTDIR"
}
export -f run_one
export SCRIPT_DIR DOMAINS_TSV TBL_DIR TIDK_DIR OUTDIR

cut -f1,2 "$ASSEMBLY_LIST" \
  | parallel --colsep '\t' -j "$JOBS" run_one {1} {2}

echo "[run_extract_tr_domains] done. Outputs: $OUTDIR"
