#!/usr/bin/env bash
set -euo pipefail

# ./src/pipeline_status.bash
#
# Reports, for every species in inputs/dtol_plant_paths.txt, which pipeline
# stages have output already and which don't - so after
# update_assembly_list.bash adds new/changed species, it's obvious what
# still needs running (and each stage's run_*.bash script already skips
# species with existing output, so the counts here are exactly the queue
# size for each `bsub`/local run).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
total=$(wc -l < "$ASSEMBLY_LIST")

count_done () {
  local suffix="$1" dir="$2"
  local n=0
  while IFS=$'\t' read -r species _; do
    [[ -z "$species" ]] && continue
    # -e not -s: TERT tblastn's output can legitimately be a completed,
    # valid, empty (0-byte) file for a genuine zero-hit species (e.g.
    # Viscum_album) - see worker_tblastn_tert.bash's atomic write.
    [[ -e "${dir}/${species}${suffix}" ]] && n=$((n + 1))
  done < "$ASSEMBLY_LIST"
  echo "$n"
}

echo "[pipeline_status] $total species in $ASSEMBLY_LIST"
echo
printf "%-38s %8s %8s\n" "stage" "done" "missing"
for stage in \
  "TR nhmmer|_all_tr_seqs_020226.tbl|${OUTPUT_DIR}/tr_tbls" \
  "TERT tblastn|.tbl|${OUTPUT_DIR}/tert_tblastn" \
  "TIDK|.tidk.tsv|${OUTPUT_DIR}/tidk" \
  "TR domain extraction|.tsv|${OUTPUT_DIR}/tr_domains" \
  "TR multi-locus extraction|.tsv|${OUTPUT_DIR}/tr_template_loci" \
  ; do
  IFS='|' read -r name suffix dir <<< "$stage"
  done_n=$(count_done "$suffix" "$dir")
  printf "%-38s %8s %8s\n" "$name" "$done_n" "$((total - done_n))"
done

echo
echo "Note: TR multi-locus extraction only applies to species with >=2 independent"
echo "TR-homologous loci, so its 'missing' count includes single-locus species by design -"
echo "see src/run_extract_tr_template_all_loci.bash, which skips those cheaply."
