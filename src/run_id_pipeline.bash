#!/usr/bin/env bash
set -euo pipefail

# ./src/run_id_pipeline.bash submit|finish|status
#
# Single wrapper for steps 3-11 of the README's "Onboarding new genomes"
# workflow (after you've already done steps 1-2: reviewed
# update_assembly_list.bash's report and run it with --apply). Every
# stage script it calls already skips species with existing output, so
# all three modes are safe to re-run at any time.
#
# Deliberately NOT a single linear "run everything and wait" script:
# run_tr_hmms.bash, run_tert_hmms.bash, run_tblastn_tert_bsub.bash, and
# run_tidk.bash all submit LSF jobs internally and return as soon as
# submission is done, not when the jobs finish - one real case of
# run_tblastn_tert_bsub.bash took 12+ days for a single large genome.
# The local extraction stages (get_all_tr_seqs.bash, recreate_TR_hmm.bash,
# run_extract_tr_domains.bash, run_extract_tr_template_all_loci.bash) need
# TR nhmmer + TIDK to have actually COMPLETED, not just been submitted -
# chaining everything into one blocking script would either hang for days
# or silently extract from incomplete data. So:
#
#   submit  - fires off the 4 independent LSF-backed search stages (TR
#             nhmmer, TERT nhmmer, TERT tblastn, TIDK) for every species
#             that doesn't have output yet. Returns immediately; the
#             actual jobs run on LSF in the background. TERT completion
#             does not gate anything else in this pipeline (TR/TIDK
#             extraction doesn't depend on it), so it's fine for TERT
#             tblastn to still be running when you move on.
#   finish  - shows current per-stage completion (via pipeline_status.bash)
#             then runs the local extraction stages in dependency order
#             (get_all_tr_seqs -> recreate_TR_hmm -> run_extract_tr_domains
#             + run_extract_tr_template_all_loci). Safe to run before TR
#             nhmmer/TIDK are 100% done - it'll just extract whatever's
#             ready and skip the rest; re-run later to pick up the rest.
#   status  - just runs pipeline_status.bash (no action) - check before
#             deciding whether to run submit again or move to finish.
#
# NOT included: map_tr_domains.py. That rebuilds the domain-boundary map
# from the TR HMM's match-state columns, only needed if the HMM profile
# itself changes (recreate_TR_hmm.bash's FORCE=1 rebuild) - a deliberate,
# rare, manual decision, not something to trigger automatically for every
# routine batch of new species.
#
# Typical flow for a batch of newly-added species:
#   bash src/run_id_pipeline.bash submit     # kick off TR/TERT/TIDK searches
#   ...come back after LSF jobs have had time to run...
#   bash src/run_id_pipeline.bash status     # check what's done
#   bash src/run_id_pipeline.bash finish     # extract domains/templates for what's ready
#   ...repeat status/finish once remaining jobs complete...

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MODE="${1:-}"
if [[ "$MODE" != "submit" && "$MODE" != "finish" && "$MODE" != "status" ]]; then
  echo "Usage: $0 submit|finish|status" >&2
  exit 1
fi

if [[ "$MODE" == "status" ]]; then
  bash "${SCRIPT_DIR}/pipeline_status.bash"
  exit 0
fi

if [[ "$MODE" == "submit" ]]; then
  echo "[run_id_pipeline] submitting TR nhmmer jobs..."
  bash "${SCRIPT_DIR}/run_tr_hmms.bash"
  echo "[run_id_pipeline] submitting TERT nhmmer jobs..."
  bash "${SCRIPT_DIR}/run_tert_hmms.bash"
  echo "[run_id_pipeline] submitting TERT tblastn jobs (can take days per genome - not a gate for anything else)..."
  bash "${SCRIPT_DIR}/run_tblastn_tert_bsub.bash"
  echo "[run_id_pipeline] submitting TIDK jobs..."
  bash "${SCRIPT_DIR}/run_tidk.bash"
  echo
  echo "[run_id_pipeline] All submission stages fired. Jobs run on LSF in the"
  echo "[run_id_pipeline] background - check progress with 'bjobs' or"
  echo "[run_id_pipeline] 'bash src/run_id_pipeline.bash status', then run"
  echo "[run_id_pipeline] 'bash src/run_id_pipeline.bash finish' once TR nhmmer"
  echo "[run_id_pipeline] and TIDK are done for the species you care about."
  exit 0
fi

# finish
echo "[run_id_pipeline] current status:"
bash "${SCRIPT_DIR}/pipeline_status.bash"
echo
echo "[run_id_pipeline] running local extraction stages (skips species whose"
echo "[run_id_pipeline] prerequisites - TR nhmmer, TIDK - aren't done yet;"
echo "[run_id_pipeline] re-run 'finish' again later to pick those up)..."
echo
echo "[run_id_pipeline] consolidating TR sequences..."
bash "${SCRIPT_DIR}/get_all_tr_seqs.bash"
echo "[run_id_pipeline] TR HMM (skips rebuild unless FORCE=1 - see script)..."
bash "${SCRIPT_DIR}/recreate_TR_hmm.bash"
echo "[run_id_pipeline] extracting TR domains (needs TR nhmmer + TIDK)..."
bash "${SCRIPT_DIR}/run_extract_tr_domains.bash"
echo "[run_id_pipeline] extracting TR templates across all loci (needs TR nhmmer + TR HMM)..."
bash "${SCRIPT_DIR}/run_extract_tr_template_all_loci.bash"
echo
echo "[run_id_pipeline] finish complete. Status now:"
bash "${SCRIPT_DIR}/pipeline_status.bash"
