#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

PATH_TO_ESLSFETCH="/software/team301/hmmer-3.4/easel/miniapps/esl-sfetch"

TR_TBL_DIR="${OUTPUT_DIR}/tr_tbls"
OUT_FASTA="${TR_TBL_DIR}/all_tr_seqs.fasta"

mkdir -p "$TR_TBL_DIR"
FORCE=1

if [[ -s "$OUT_FASTA" && "${FORCE:-0}" != "1" ]]; then
  echo "[get_all_tr_seqs] $OUT_FASTA exists and is non-empty; skipping. Set FORCE=1 to overwrite."
  exit 0
fi

: > "$OUT_FASTA"

for i in "${TR_TBL_DIR}"/*.tbl; do
  [[ ! -f "$i" ]] && continue
  BN=$(basename "$i" .tbl)
  echo "[get_all_tr_seqs] Extracting from $i"
  extract_nhmmer_tblout -v 0.00001 -e "$PATH_TO_ESLSFETCH" -s "$BN" "$i" >> "$OUT_FASTA"
done

echo "[get_all_tr_seqs] Combined TR sequences written to $OUT_FASTA"

