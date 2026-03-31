#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
WORKER="${WORKER:-${SCRIPT_DIR}/worker_tblastn_tert.bash}"

export BLASTBIN="${BLASTBIN:-/software/team301/ncbi-blast-2.16.0+/bin}"
export THREADS="${THREADS:-8}"
export QUERY_AA="${QUERY_AA:-${INPUT_DIR}/plant_TERT_regions.faa}"

export OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tert_tblastn}"
export LOGDIR="${LOGDIR:-${OUTPUT_DIR}/logs/tert_tblastn}"
export EVALUE="${EVALUE:-1e-4}"
export MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-10}"

mkdir -p "$OUTDIR" "$LOGDIR"

BSUB_THREADS="${BSUB_THREADS:-8}"
BSUB_MEM_GB="${BSUB_MEM_GB:-60}"
BSUB_QUEUE="${BSUB_QUEUE:-basement}"

BSUB_MEM_MB=$(( BSUB_MEM_GB * 1024 ))

submit_one () {
  local species="$1" fasta="$2"
  local out="${OUTDIR}/${species}.tbl"
  local jobname="tblastn_${species}"

  if [[ -s "$out" ]]; then
    echo "[driver-tblastn] ${species}: $out exists, skipping."
    return
  fi

  export THREADS="$BSUB_THREADS"

  local log_out="${LOGDIR}/${species}.%J.out"
  local log_err="${LOGDIR}/${species}.%J.err"
  local queue_flag=()
  [[ -n "$BSUB_QUEUE" ]] && queue_flag=(-q "$BSUB_QUEUE")

  bsub \
    -J "$jobname" \
    "${queue_flag[@]}" \
    -n "$BSUB_THREADS" \
    -R "span[hosts=1] select[mem>=${BSUB_MEM_MB}] rusage[mem=${BSUB_MEM_MB}]" \
    -M "$BSUB_MEM_MB" \
    -o "$log_out" \
    -e "$log_err" \
    -env "all" \
    bash -lc "\"$WORKER\" \"$species\" \"$fasta\""
}

cut -f1,2 "$ASSEMBLY_LIST" | while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue
  submit_one "$species" "$fasta"
done

echo "[run_tblastn_tert_bsub] Submitted tBLASTn TERT jobs. Outputs: $OUTDIR  Logs: $LOGDIR"

