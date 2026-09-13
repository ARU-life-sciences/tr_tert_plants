#!/usr/bin/env bash
set -uo pipefail

# ./src/run_confirm_differing_sets.bash [core_template_sets.tsv]
#
# For every "differing" comparable set from tr_core_template_sets.py (i.e.
# a species where >1 distinct core-template sequence exists among its TR
# loci), runs confirm_repeat_positional.py for BOTH the dominant and the
# top minority variant - the question being: does the specific
# point-mutant sequence carried by the minority TR gene copies actually
# show up at the species' chromosome termini, the same way the dominant
# variant does?
#
# Deduplicates (species, variant) pairs across sets before running, since
# a species can appear in several differing sets (different core
# lengths) and dominant/minority strings can repeat. Local (no bsub) -
# each confirm_repeat_positional.py call runs `tidk search` directly on
# the (possibly gzipped) genome, so no separate decompression step is
# needed here. See run_confirm_differing_sets_bsub.bash for the LSF
# version when running the full dataset.
#
# Safe to re-run: skips (species, variant) pairs with existing output.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

SETS_TSV="${1:-${OUTPUT_DIR}/tr_repeat_correlation/core_template_sets.tsv}"
OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tr_repeat_correlation/positional}"
JOBS="${JOBS:-10}"
mkdir -p "$OUTDIR"

run_variant () {
  species="$1"; variant="$2"
  outfile="${OUTDIR}/${species}_${variant}_telomeric_repeat_windows.tsv"
  if [[ -s "$outfile" ]]; then
    echo "[run_confirm_differing_sets] ${species} ${variant}: exists, skipping."
    return
  fi
  echo "[run_confirm_differing_sets] ${species} ${variant}: running"
  python3 "${SCRIPT_DIR}/confirm_repeat_positional.py" "$species" "$variant" "$OUTDIR" > /dev/null
}
export -f run_variant
export OUTDIR SCRIPT_DIR

# columns: species core_length n_loci n_distinct dominant dom_n dom_target minority min_n min_target
awk -F'\t' 'NR>1 && $4>1 {print $1"\t"$5; print $1"\t"$8}' "$SETS_TSV" | sort -u \
  | parallel --colsep '\t' -j "$JOBS" run_variant {1} {2}

echo "[run_confirm_differing_sets] done. Outputs: $OUTDIR"
