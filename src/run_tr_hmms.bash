#!/usr/bin/env bash
set -euo pipefail

# ./src/run_tr_hmms.bash
# note this can be re-run after we have done an initial round
# of identification. In our case, the new and updated hmm file
# for the TR gene is in:
# ../outputs/tr_tbls/all_tr_seqs.hmm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

TR_HMM="${INPUT_DIR}/TR.hmm"
UPDATED_TR_HMM="${OUTPUT_DIR}/tr_tbls/all_tr_seqs.hmm"
ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"

OUTDIR="${OUTPUT_DIR}/tr_tbls"
LOGDIR="${OUTPUT_DIR}/logs/TR"
mkdir -p "$OUTDIR" "$LOGDIR"

mbMem=10000    # 10 GB

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue

  out_tbl="${OUTDIR}/${species}_all_tr_seqs_020226.tbl"
  if [[ -s "$out_tbl" ]]; then
    echo "[TR nhmmer] ${species}: $out_tbl exists, skipping."
    continue
  fi

  echo "[TR nhmmer] Submitting $species"
  bsub -n 20 -q normal \
    -R"span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M"${mbMem}" \
    -o "${LOGDIR}/${species}.out" \
    -e "${LOGDIR}/${species}.err" \
    "/software/team301/hmmer-3.4/src/nhmmer --dna --cpu 20 \
       --tblout \"${out_tbl}\" \
       \"${UPDATED_TR_HMM}\" \"${fasta}\""

done < "$ASSEMBLY_LIST"

echo "[run_tr_hmms] Submitted TR nhmmer jobs. Outputs: $OUTDIR"

