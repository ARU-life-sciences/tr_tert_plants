#!/usr/bin/env bash
set -euo pipefail

# ./src/run_extract_tr_template_all_loci_bsub.bash
#
# LSF counterpart to run_extract_tr_template_all_loci.bash - submits one
# bsub job per multi-locus species instead of running locally via GNU
# parallel. Depends on that species' TR nhmmer job ("tr_<species>") via
# -w when the tbl doesn't exist yet, same pattern as
# run_extract_tr_domains_bsub.bash.
#
# The multi-locus check itself (>=2 decent loci) needs the tbl to already
# exist, so for genuinely new species (no tbl yet) this submits
# unconditionally whenever a tr_<species> job was just submitted, and lets
# the job itself exit quietly if it turns out there's only one locus.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
TBL_DIR="${OUTPUT_DIR}/tr_tbls"
HMM="${HMM:-${TBL_DIR}/all_tr_seqs.hmm}"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_template_loci}"
LOGDIR="${OUTPUT_DIR}/logs/tr_template_loci"
mkdir -p "$OUTDIR" "$LOGDIR"

mbMem=4000  # 4 GB

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue

  out="${OUTDIR}/${species}.tsv"
  if [[ -s "$out" ]]; then
    echo "[run_extract_tr_template_all_loci_bsub] ${species}: $out exists, skipping."
    continue
  fi

  tbl="${TBL_DIR}/${species}_all_tr_seqs_020226.tbl"

  if [[ -s "$tbl" ]]; then
    # tbl already exists - do the cheap >=2-loci pre-check now rather than
    # submitting a job for the (majority) single-locus species.
    n=$(awk '!/^#/ && $13<=1e-4 && ($6-$5)>=100' "$tbl" | wc -l)
    if [[ "$n" -lt 2 ]]; then
      echo "[run_extract_tr_template_all_loci_bsub] ${species}: only $n decent loci, skipping."
      continue
    fi
    wait_clause=()
  else
    # tbl doesn't exist yet - a tr_<species> job was (presumably) just
    # submitted by run_tr_hmms.bash; depend on it and let the job itself
    # decide whether there are enough loci once the tbl exists.
    wait_clause=(-w "done(tr_${species})")
  fi

  echo "[run_extract_tr_template_all_loci_bsub] Submitting ${species} ${wait_clause[*]:+(waiting on: ${wait_clause[1]})}"
  bsub -J "trloci_${species}" "${wait_clause[@]}" -n 2 -q normal \
    -R "span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M "${mbMem}" \
    -o "${LOGDIR}/${species}.out" \
    -e "${LOGDIR}/${species}.err" \
    "python3 '${SCRIPT_DIR}/extract_tr_template_all_loci.py' '${species}' '${fasta}' '${tbl}' '${HMM}' '${OUTDIR}'"

done < "$ASSEMBLY_LIST"

echo "[run_extract_tr_template_all_loci_bsub] Submitted. Outputs: $OUTDIR"
