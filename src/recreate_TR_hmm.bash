#!/usr/bin/env bash
set -euo pipefail
# using the `../working/TR_tbl_outputs/all_tr_seqs.fasta`
# recreate the HMM, which was initially formed from fajkus.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

TR_TBL_DIR="${OUTPUT_DIR}/tr_tbls"

ALN="${TR_TBL_DIR}/all_tr_seqs.aln"
HMM="${TR_TBL_DIR}/all_tr_seqs.hmm"

if [[ -s "$HMM" && "${FORCE:-0}" != "1" ]]; then
  echo "[recreate_TR_hmm] $HMM already exists; skipping. Set FORCE=1 to rebuild."
  exit 0
fi

# ALIGN (uncomment / tweak MAFFT settings as needed)
/software/team301/mafft-7.525-with-extensions/core/mafft \
  --genafpair --maxiterate 1000 \
  "${TR_TBL_DIR}/all_tr_seqs.fasta" > "$ALN"

# BUILD HMM
/software/team301/hmmer-3.4/src/hmmbuild -n all_tr_seqs "$HMM" "$ALN"
software/team301/hmmer-3.4/src/hmmpress "$HMM"

echo "[recreate_TR_hmm] HMM written to $HMM"

