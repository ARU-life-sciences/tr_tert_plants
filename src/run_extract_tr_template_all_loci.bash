#!/usr/bin/env bash
set -euo pipefail

# ./src/run_extract_tr_template_all_loci.bash
#
# Runs extract_tr_template_all_loci.py for every species with >=2 decent
# (E<=1e-4, hmm span>=100 cols) TR-homologous loci in
# outputs/tr_tbls/<species>_all_tr_seqs_020226.tbl, to assess within-species
# Template variation across paralogous/duplicated loci (motivated by the
# Carlina vulgaris finding - see notes/carlina_vulgaris_telomere_repeat.md).
#
# Local (no bsub), parallelised like run_extract_tr_domains.bash. Safe to
# re-run: skips species that already have an output file.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
TBL_DIR="${OUTPUT_DIR}/tr_tbls"
HMM="${HMM:-${TBL_DIR}/all_tr_seqs.hmm}"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_template_loci}"
JOBS="${JOBS:-12}"

mkdir -p "$OUTDIR"

run_one () {
  species="$1"; fasta="$2"
  out="${OUTDIR}/${species}.tsv"
  if [[ -s "$out" ]]; then
    echo "[run_extract_tr_template_all_loci] ${species}: $out exists, skipping."
    return
  fi
  tbl="${TBL_DIR}/${species}_all_tr_seqs_020226.tbl"
  if [[ ! -s "$tbl" ]]; then
    echo "[run_extract_tr_template_all_loci] ${species}: no tbl at $tbl, skipping."
    return
  fi
  # Cheap pre-check: skip the (majority) single-locus species without
  # touching the genome at all.
  n=$(awk '!/^#/ && $13<=1e-4 && ($6-$5)>=100' "$tbl" | wc -l)
  if [[ "$n" -lt 2 ]]; then
    echo "[run_extract_tr_template_all_loci] ${species}: only $n decent loci, skipping."
    return
  fi
  python3 "${SCRIPT_DIR}/extract_tr_template_all_loci.py" "$species" "$fasta" "$tbl" "$HMM" "$OUTDIR"
}
export -f run_one
export SCRIPT_DIR HMM TBL_DIR OUTDIR

cut -f1,2 "$ASSEMBLY_LIST" \
  | parallel --colsep '\t' -j "$JOBS" run_one {1} {2}

echo "[run_extract_tr_template_all_loci] done. Outputs: $OUTDIR"
