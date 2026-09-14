#!/usr/bin/env bash
set -euo pipefail

# Weekly cron entry point - runs update_assembly_list.bash in report-only
# mode (never --apply: this must never silently change
# inputs/dtol_plant_paths.txt or touch any pipeline outputs unattended)
# and appends a timestamped block to outputs/logs/update_assembly_list_cron.log.
#
# Installed via crontab -l; see notes/ for the onboarding workflow this
# feeds into (README.md "Onboarding new genomes"). This script only
# reports what changed - a human still runs
# `bash src/update_assembly_list.bash --apply` (and the relevant
# pipeline stages) after reviewing the log.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG="${PROJECT_ROOT}/outputs/logs/update_assembly_list_cron.log"

{
  echo "===== $(date -Iseconds) ====="
  bash "${SCRIPT_DIR}/update_assembly_list.bash" || echo "[cron_check_new_species] update_assembly_list.bash exited non-zero"
  echo
} >> "$LOG" 2>&1
