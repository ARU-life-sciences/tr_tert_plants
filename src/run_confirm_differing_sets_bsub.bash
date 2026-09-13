#!/usr/bin/env bash
set -uo pipefail

# ./src/run_confirm_differing_sets_bsub.bash [core_template_sets.tsv]
#
# LSF counterpart to run_confirm_differing_sets.bash - submits one bsub
# job per (species, variant) pair instead of running locally via GNU
# parallel. Same dedup logic and skip-if-exists behaviour.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

SETS_TSV="${1:-${OUTPUT_DIR}/tr_repeat_correlation/core_template_sets.tsv}"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_repeat_correlation/positional}"
LOGDIR="${OUTPUT_DIR}/logs/tr_repeat_correlation_positional"
mkdir -p "$OUTDIR" "$LOGDIR"

mbMem=8000  # 8 GB - tidk search reads the whole (possibly gzipped) genome

awk -F'\t' 'NR>1 && $4>1 {print $1"\t"$5; print $1"\t"$8}' "$SETS_TSV" | sort -u | \
while IFS=$'\t' read -r species variant; do
  [[ -z "$species" || -z "$variant" ]] && continue

  outfile="${OUTDIR}/${species}_${variant}_telomeric_repeat_windows.tsv"
  if [[ -s "$outfile" ]]; then
    echo "[run_confirm_differing_sets_bsub] ${species} ${variant}: exists, skipping."
    continue
  fi

  jobname="confirmrepeat_${species}_${variant}"
  echo "[run_confirm_differing_sets_bsub] Submitting ${species} ${variant}"
  bsub -J "$jobname" -n 2 -q normal \
    -R "span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M "${mbMem}" \
    -o "${LOGDIR}/${species}_${variant}.out" \
    -e "${LOGDIR}/${species}_${variant}.err" \
    "python3 '${SCRIPT_DIR}/confirm_repeat_positional.py' '${species}' '${variant}' '${OUTDIR}'"

done

echo "[run_confirm_differing_sets_bsub] Submitted. Outputs: $OUTDIR"
