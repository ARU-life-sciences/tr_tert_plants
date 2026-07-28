#!/usr/bin/env bash
set -euo pipefail

# Non-LSF counterpart to run_tblastn_tert_bsub.bash: runs worker_tblastn_tert.bash
# directly (in series, on the local machine) instead of submitting bsub jobs.
# Intended for small assembly lists that don't warrant the cluster -- e.g.
# inputs/extra_genomes.txt, populated by fetch_extra_genomes.bash -- rather
# than the full DToL run.
#
# Output lands in the same place as the bsub path (outputs/tert_tblastn/<species>.tbl),
# so extra genomes show up alongside the DToL species in any downstream analysis.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/extra_genomes.txt}"
WORKER="${WORKER:-${SCRIPT_DIR}/worker_tblastn_tert.bash}"

export QUERY_AA="${QUERY_AA:-${INPUT_DIR}/plant_TERT_regions.faa}"
export THREADS="${THREADS:-4}"

# worker_tblastn_tert.bash defaults BLASTBIN to the Sanger HPC install; fall
# back to whatever tblastn is on PATH when that doesn't exist (e.g. local runs).
if [[ -z "${BLASTBIN:-}" ]] && [[ ! -x "/software/team301/ncbi-blast-2.16.0+/bin/tblastn" ]]; then
  export BLASTBIN="$(dirname "$(command -v tblastn)")"
fi

cut -f1,2 "$ASSEMBLY_LIST" | while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue
  echo "[run_tblastn_tert_local] ${species}"
  "$WORKER" "$species" "$fasta"
done

echo "[run_tblastn_tert_local] Done."
