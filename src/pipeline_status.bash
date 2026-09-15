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

# TERT HMM search (run_tert_hmms.bash) writes 3 per-species files in a
# per-species subdirectory, not one file matching a simple suffix - needs
# its own check, matching that script's own skip-condition exactly.
count_tert_hmm_done () {
  local dir="${OUTPUT_DIR}/tert_tbls"
  local n=0
  while IFS=$'\t' read -r species _; do
    [[ -z "$species" ]] && continue
    local spp_dir="${dir}/${species}"
    if [[ -s "${spp_dir}/${species}_isoformx1.tbl" && \
          -s "${spp_dir}/${species}_catalytic_subunit.tbl" && \
          -s "${spp_dir}/${species}_tert1.tbl" ]]; then
      n=$((n + 1))
    fi
  done < "$ASSEMBLY_LIST"
  echo "$n"
}

print_row () {
  local name="$1" done_n="$2"
  printf "%-38s %8s %8s\n" "$name" "$done_n" "$((total - done_n))"
}

echo "[pipeline_status] $total species in $ASSEMBLY_LIST"
echo
printf "%-38s %8s %8s\n" "stage" "done" "missing"
print_row "TR nhmmer" "$(count_done "_all_tr_seqs_020226.tbl" "${OUTPUT_DIR}/tr_tbls")"
print_row "TERT HMM (nhmmer)" "$(count_tert_hmm_done)"
print_row "TERT tblastn" "$(count_done ".tbl" "${OUTPUT_DIR}/tert_tblastn")"
print_row "TIDK" "$(count_done ".tidk.tsv" "${OUTPUT_DIR}/tidk")"
print_row "TR domain extraction" "$(count_done ".tsv" "${OUTPUT_DIR}/tr_domains")"
print_row "TR multi-locus extraction" "$(count_done ".tsv" "${OUTPUT_DIR}/tr_template_loci")"

echo
echo "Note: TR multi-locus extraction only applies to species with >=2 independent"
echo "TR-homologous loci, so its 'missing' count includes single-locus species by design -"
echo "see src/run_extract_tr_template_all_loci.bash, which skips those cheaply."
