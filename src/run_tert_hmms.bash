#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tert_hmms.bash
#
# Submit nhmmer jobs for TERT HMMs across all plant assemblies using LSF.
# - Reads a list of assemblies (species<TAB>path) from:
#       inputs/dtol_plant_paths.txt     (overridable with $ASSEMBLY_LIST)
# - Uses HMMs in:
#       inputs/clusters/*.hmm
# - Writes outputs to:
#       outputs/tert_tbls/<species>/*.tbl
# - Skips species where ALL expected .tbl files already exist and are non-empty.
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

# Assembly list: species<TAB>/path/to/assembly.fa[.gz]
ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"

# HMMs
HMM_DIR="${INPUT_DIR}/clusters"
ISOFORM_X1_HMM="${HMM_DIR}/isoform_X1_tert.hmm"
CATALYTIC_HMM="${HMM_DIR}/catalytic_subunit_tert.hmm"
TERT1_HMM="${HMM_DIR}/tert1.hmm"

# Outputs + logs
OUTBASE="${OUTPUT_DIR}/tert_tbls"
LOGDIR="${OUTPUT_DIR}/logs/tert_nhmmer"
mkdir -p "$OUTBASE" "$LOGDIR"

# LSF resources
mbMem="${mbMem:-20000}"         # 20 GB default; override with mbMem=...
nCPU="${nCPU:-20}"
QUEUE="${QUEUE:-small}"

echo "[run_tert_hmms] Using assembly list: $ASSEMBLY_LIST"
echo "[run_tert_hmms] HMMs: "
echo "  $ISOFORM_X1_HMM"
echo "  $CATALYTIC_HMM"
echo "  $TERT1_HMM"
echo "[run_tert_hmms] Output base: $OUTBASE"
echo "[run_tert_hmms] Logs: $LOGDIR"
echo

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue

  spp_dir="${OUTBASE}/${species}"
  mkdir -p "$spp_dir"

  iso_tbl="${spp_dir}/${species}_isoformx1.tbl"
  cat_tbl="${spp_dir}/${species}_catalytic_subunit.tbl"
  tert1_tbl="${spp_dir}/${species}_tert1.tbl"

  # Skip if all outputs already exist and are non-empty
  if [[ -s "$iso_tbl" && -s "$cat_tbl" && -s "$tert1_tbl" ]]; then
    echo "[run_tert_hmms] ${species}: all TERT tbls present, skipping."
    continue
  fi

  log_out="${LOGDIR}/${species}.%J.out"
  log_err="${LOGDIR}/${species}.%J.err"
  jobname="tert_nhmmer_${species}"

  echo "[run_tert_hmms] Submitting ${species}"
  echo "  assembly: $fasta"
  echo "  outputs :"
  echo "    $iso_tbl"
  echo "    $cat_tbl"
  echo "    $tert1_tbl"
  echo

  # One LSF job per species, running all three nhmmer searches in series
  bsub \
    -J "$jobname" \
    -q "$QUEUE" \
    -n "$nCPU" \
    -R"span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M"$mbMem" \
    -o "$log_out" \
    -e "$log_err" \
    "/software/team301/hmmer-3.4/src/nhmmer --dna --cpu ${nCPU} \
        --tblout \"${iso_tbl}\" \"${ISOFORM_X1_HMM}\" \"${fasta}\" && \
     /software/team301/hmmer-3.4/src/nhmmer --dna --cpu ${nCPU} \
        --tblout \"${cat_tbl}\" \"${CATALYTIC_HMM}\" \"${fasta}\" && \
     /software/team301/hmmer-3.4/src/nhmmer --dna --cpu ${nCPU} \
        --tblout \"${tert1_tbl}\" \"${TERT1_HMM}\" \"${fasta}\""

done < "$ASSEMBLY_LIST"

echo "[run_tert_hmms] Submitted TERT nhmmer jobs."

