#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

species="$1"
fasta_path="$2"

BLASTBIN="${BLASTBIN:-/software/team301/ncbi-blast-2.16.0+/bin}"
THREADS="${THREADS:-8}"
QUERY_AA="${QUERY_AA:-${PROJECT_ROOT}/inputs/plant_TERT_regions.faa}"

OUTDIR="${OUTDIR:-${OUTPUT_DIR}/tert_tblastn}"
LOGDIR="${LOGDIR:-${OUTPUT_DIR}/logs/tert_tblastn}"
EVALUE="${EVALUE:-1e-4}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-10}"

mkdir -p "$OUTDIR" "$LOGDIR"

out="${OUTDIR}/${species}.tbl"
if [[ -s "$out" ]]; then
  echo "[worker-tblastn] ${species}: $out exists, skipping."
  exit 0
fi

tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}" "tblastn_${species}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

fa="${tmpdir}/${species}.fa"
dbpref="${tmpdir}/${species}"

echo "[worker-tblastn] ${species}: preparing FASTA"
if [[ "$fasta_path" =~ \.gz$ ]]; then
  gunzip -c "$fasta_path" > "$fa"
else
  rsync -a "$fasta_path" "$fa"
fi

echo "[worker-tblastn] ${species}: makeblastdb"
"$BLASTBIN/makeblastdb" -in "$fa" -dbtype nucl -parse_seqids -out "$dbpref" \
  1>"${LOGDIR}/${species}.tblastn.makeblastdb.log" 2>&1

echo "[worker-tblastn] ${species}: tblastn"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
"$BLASTBIN/tblastn" \
  -query "$QUERY_AA" \
  -db "$dbpref" \
  -num_threads "$THREADS" \
  -seg yes \
  -evalue "$EVALUE" \
  -max_target_seqs "$MAX_TARGET_SEQS" \
  -comp_based_stats 1 \
  -outfmt "$FMT" \
  -out "$out"

echo "[worker-tblastn] ${species}: done."

