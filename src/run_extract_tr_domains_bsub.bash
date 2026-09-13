#!/usr/bin/env bash
set -euo pipefail

# ./src/run_extract_tr_domains_bsub.bash
#
# LSF counterpart to run_extract_tr_domains.bash (which runs locally via
# GNU parallel) - submits one bsub job per species instead. Each job
# depends on that species' TR nhmmer job ("tr_<species>", from
# run_tr_hmms.bash) and TIDK job ("tidk_<species>", from run_tidk.bash)
# via LSF's -w, so it's safe to submit this in the same breath as those
# two even though it needs their output: LSF just holds it PEND until both
# report done. If a species' tbl/tidk output already exists from an
# earlier run (no new job submitted this time), -w is simply omitted and
# the job runs immediately.
#
# Skips species with existing output, same as the local version.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
TBL_DIR="${OUTPUT_DIR}/tr_tbls"
DOMAINS_TSV="${DOMAINS_TSV:-${INPUT_DIR}/TR_domains.all_tr_seqs.tsv}"
TIDK_DIR="${OUTPUT_DIR}/tidk"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_domains}"
LOGDIR="${OUTPUT_DIR}/logs/tr_domains"
mkdir -p "$OUTDIR" "$LOGDIR"

mbMem=4000  # 4 GB - dominated by one genome decompression + a small hmmalign

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue

  out="${OUTDIR}/${species}.tsv"
  if [[ -s "$out" ]]; then
    echo "[run_extract_tr_domains_bsub] ${species}: $out exists, skipping."
    continue
  fi

  tbl="${TBL_DIR}/${species}_all_tr_seqs_020226.tbl"
  tidk="${TIDK_DIR}/${species}.tidk.tsv"
  [[ -s "$tidk" ]] || tidk=""

  # Only add a -w dependency for stages that were actually just submitted
  # (i.e. their output doesn't exist yet) - bjobs -J matches on jobs LSF
  # still knows about, but it's simplest and robust to just always depend
  # on the job name when the tbl/tidk file is missing, and skip -w when it
  # already exists (nothing to wait for).
  wait_clause=()
  conditions=()
  [[ ! -s "$tbl" ]] && conditions+=("done(tr_${species})")
  [[ -z "$tidk" ]] && conditions+=("done(tidk_${species})")
  if [[ ${#conditions[@]} -gt 0 ]]; then
    joined="${conditions[0]}"
    for c in "${conditions[@]:1}"; do
      joined="${joined} && ${c}"
    done
    wait_clause=(-w "$joined")
  fi

  echo "[run_extract_tr_domains_bsub] Submitting ${species} ${wait_clause[*]:+(waiting on: ${wait_clause[1]})}"
  bsub -J "trdomains_${species}" "${wait_clause[@]}" -n 2 -q normal \
    -R "span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M "${mbMem}" \
    -o "${LOGDIR}/${species}.out" \
    -e "${LOGDIR}/${species}.err" \
    "python3 '${SCRIPT_DIR}/extract_tr_domains.py' '${species}' '${fasta}' '${tbl}' '${DOMAINS_TSV}' '${tidk}' '${OUTDIR}'"

done < "$ASSEMBLY_LIST"

echo "[run_extract_tr_domains_bsub] Submitted. Outputs: $OUTDIR"
