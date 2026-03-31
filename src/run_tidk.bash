#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tidk_explore.bash
#
# Run "tidk explore" across all plant assemblies using LSF.
# Reads species<TAB>path from: inputs/dtol_plant_paths.txt
# Outputs: outputs/tidk/<species>.tidk.tsv
# Skips species with an existing output (unless FORCE=1).
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

# List of species + genome paths
ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"

# Output dir
OUTDIR="${OUTPUT_DIR}/tidk"
LOGDIR="${OUTPUT_DIR}/logs/tidk"
mkdir -p "$OUTDIR" "$LOGDIR"

# TIDK binary
TIDK="${TIDK:-tidk}"

# Default parameters for tidk explore
MINLEN="${MINLEN:-4}"      # -m
MAXLEN="${MAXLEN:-30}"     # -x
THRESH="${THRESH:-100}"    # --threshold (tidk default)
DIST="${DIST:-0.01}"       # --distance

# LSF resources (override as needed)
nCPU="${nCPU:-8}"
memGB="${memGB:-8}"
QUEUE="${QUEUE:-small}"

memMB=$(( memGB * 1024 ))

echo "[run_tidk_explore] assembly list: $ASSEMBLY_LIST"
echo "[run_tidk_explore] output dir:   $OUTDIR"
echo "[run_tidk_explore] logs:         $LOGDIR"
echo "[run_tidk_explore] tidk params:  -m $MINLEN -x $MAXLEN --threshold $THRESH --distance $DIST"
echo

while IFS=$'\t' read -r species fasta; do
    [[ -z "$species" || -z "$fasta" ]] && continue

    out="${OUTDIR}/${species}.tidk.tsv"

    # Skip if output already present unless FORCE=1
    if [[ -s "$out" && "${FORCE:-0}" != "1" ]]; then
        echo "[tidk] ${species}: output exists, skipping."
        continue
    fi

    echo "[tidk] Submitting: ${species}"
    jobname="tidk_${species}"
    log_out="${LOGDIR}/${species}.%J.out"
    log_err="${LOGDIR}/${species}.%J.err"

    # One job per species
    bsub \
      -J "$jobname" \
      -q "$QUEUE" \
      -n "$nCPU" \
      -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
      -M "$memMB" \
      -o "$log_out" \
      -e "$log_err" \
      bash -lc "
        ${TIDK} explore \
            -m ${MINLEN} -x ${MAXLEN} \
            --threshold ${THRESH} \
            --distance ${DIST} \
            --log \
            '${fasta}' \
            > '${out}'
      "
done < "$ASSEMBLY_LIST"

echo "[run_tidk_explore] Submitted all jobs."

