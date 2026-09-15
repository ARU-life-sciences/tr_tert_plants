#!/usr/bin/env bash
set -euo pipefail

# ./src/update_assembly_list.bash [--apply] [--purge-stale]
#
# Refreshes inputs/dtol_plant_paths.txt against the current DToL assembly
# tree and reports what changed, rather than silently overwriting it:
#
#   NEW      species not previously in the list - safe to just add, the
#            rest of the pipeline (run_tr_hmms.bash, run_tidk.bash,
#            run_tblastn_tert_bsub.bash, run_extract_tr_domains.bash, ...)
#            already skips species with existing output, so it will only
#            do work for these.
#   CHANGED  species whose "latest primary.fa.gz" path differs from what's
#            in the list - almost always because the assembly was
#            re-curated/re-scaffolded since we last ran (this bit us twice:
#            Lycopus_europaeus and Galium_aparine both had genuinely
#            different sequence content under the same species name, not
#            just a renamed file - see git log). Downstream per-species
#            outputs are keyed by species name and skip-if-exists, so they
#            will NOT automatically pick up a changed genome; the old
#            outputs would silently keep using stale coordinates.
#   REMOVED  species that dropped out of the assembly tree entirely (e.g.
#            withdrawn). Never acted on automatically.
#
# Without --apply, this is report-only (writes the fresh list to a .new
# file for inspection, doesn't touch inputs/dtol_plant_paths.txt or any
# outputs). With --apply, updates inputs/dtol_plant_paths.txt (adds NEW,
# updates CHANGED paths). With --apply --purge-stale, additionally deletes
# the existing per-species output files for CHANGED species (tr_tbls,
# tert_tbls, tert_tblastn, tidk, tr_domains, tr_template_loci) so the next
# pipeline run regenerates them against the new genome - CHANGED species'
# outputs are otherwise left in place (stale) since we can't tell from the
# path alone whether the underlying sequence actually differs.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"
LATEST_TOOL="${LATEST_TOOL:-/software/team301/user/mb39/scripts/any_latest_assemblies.bash}"
DTOL_GROUPS=(dicots vascular-plants non-vascular-plants monocots)

APPLY=0
PURGE_STALE=0
for arg in "$@"; do
  case "$arg" in
    --apply) APPLY=1 ;;
    --purge-stale) PURGE_STALE=1 ;;
    *) echo "Unknown argument: $arg" >&2; exit 1 ;;
  esac
done
if [[ "$PURGE_STALE" == "1" && "$APPLY" != "1" ]]; then
  echo "--purge-stale requires --apply" >&2
  exit 1
fi

RAW_FRESH=$(mktemp)
bash "$LATEST_TOOL" "${DTOL_GROUPS[@]}" | sort > "$RAW_FRESH"

# Guard against organelle-only assemblies silently shadowing the real
# nuclear genome. $LATEST_TOOL picks whichever curated directory has the
# newest mtime, with no awareness of assembly type - if a species' plastid
# assembly gets curated more recently than its nuclear one (both are
# genuinely "latest" by mtime, just for different things), it wins, and TR
# nhmmer/tblastn against it just produces uniform zero-hit results that
# look exactly like a real negative finding until you check by hand. This
# bit Galanthus_nivalis (see notes/) - every downstream stage ran against
# ~157KB of plastid sequence and "found" nothing, for three whole pipeline
# runs, before anyone noticed. Two checks: an obvious naming smell, and (as
# a backstop for assemblies that don't announce it in the path) a size
# floor comfortably below any real nuclear plant genome but well above a
# compressed plastid/mitochondrial genome.
MIN_BYTES="${MIN_BYTES:-5000000}"  # 5MB
FRESH="${OUTPUT_DIR}/dtol_plant_paths.fresh.txt"
SUSPECT="${OUTPUT_DIR}/dtol_plant_paths.suspect_organelle.txt"
: > "$FRESH"
: > "$SUSPECT"
while IFS=$'\t' read -r species path; do
  [[ -z "$species" ]] && continue
  reason=""
  curdir=$(basename "$(dirname "$path")")
  if echo "$path" | grep -qiE 'plastid|chloroplast|mito(chondri)?'; then
    reason="path mentions plastid/chloroplast/mito"
  elif ! echo "$curdir" | grep -qE '^[A-Za-z]+[0-9]+(\.hap[0-9]+)?\.[0-9]+$'; then
    # Catches cobiont directories named after an organism that doesn't
    # happen to match the plastid/chloroplast/mito keyword check above
    # (e.g. laLemMinu1.Leptothrix_sp_1.1 - see notes/ for the case that
    # motivated this) - a real nuclear curated directory is always just
    # PREFIX<digits>(.hap<digits>)?.<digits>, nothing else.
    reason="curated-directory name '${curdir}' doesn't match the plain assembly-version pattern"
  elif [[ -f "$path" ]]; then
    sz=$(stat -c%s "$path" 2>/dev/null || echo 0)
    if [[ "$sz" -gt 0 && "$sz" -lt "$MIN_BYTES" ]]; then
      reason="only ${sz} bytes (< ${MIN_BYTES})"
    fi
  fi
  if [[ -n "$reason" ]]; then
    echo -e "${species}\t${path}\t${reason}" >> "$SUSPECT"
  else
    echo -e "${species}\t${path}" >> "$FRESH"
  fi
done < "$RAW_FRESH"
sort -o "$FRESH" "$FRESH"
rm -f "$RAW_FRESH"

n_suspect=$(wc -l < "$SUSPECT")
if [[ "$n_suspect" -gt 0 ]]; then
  echo "[update_assembly_list] SUSPECT (excluded, likely organelle-only - not nuclear): $n_suspect"
  sed 's/^/  ! /' "$SUSPECT"
  echo "[update_assembly_list] Not added/updated automatically. Investigate manually (see notes/) before"
  echo "[update_assembly_list] including any of these - check whether a real nuclear assembly exists yet."
  echo
fi

CUR_SORTED=$(mktemp)
sort "$ASSEMBLY_LIST" > "$CUR_SORTED"

NEW=$(comm -13 <(cut -f1 "$CUR_SORTED") <(cut -f1 "$FRESH"))
REMOVED=$(comm -23 <(cut -f1 "$CUR_SORTED") <(cut -f1 "$FRESH"))

CHANGED=""
while IFS=$'\t' read -r species path; do
  [[ -z "$species" ]] && continue
  cur_path=$(awk -F'\t' -v s="$species" '$1==s{print $2}' "$CUR_SORTED")
  if [[ -n "$cur_path" && "$cur_path" != "$path" ]]; then
    CHANGED+="${species}\n"
  fi
done < "$FRESH"
CHANGED=$(echo -e "$CHANGED" | sed '/^$/d')

n_new=$(echo -n "$NEW" | grep -c . || true)
n_changed=$(echo -n "$CHANGED" | grep -c . || true)
n_removed=$(echo -n "$REMOVED" | grep -c . || true)

echo "[update_assembly_list] current list: $(wc -l < "$CUR_SORTED") species"
echo "[update_assembly_list] fresh list:    $(wc -l < "$FRESH") species"
echo "[update_assembly_list] NEW:     $n_new"
[[ "$n_new" -gt 0 ]] && echo "$NEW" | sed 's/^/  + /'
echo "[update_assembly_list] CHANGED: $n_changed"
[[ "$n_changed" -gt 0 ]] && echo "$CHANGED" | sed 's/^/  ~ /'
echo "[update_assembly_list] REMOVED: $n_removed"
[[ "$n_removed" -gt 0 ]] && echo "$REMOVED" | sed 's/^/  - /'

if [[ "$APPLY" != "1" ]]; then
  echo
  echo "[update_assembly_list] Report-only (no --apply). Fresh list written to: $FRESH"
  echo "[update_assembly_list] Re-run with --apply to update $ASSEMBLY_LIST (NEW + CHANGED paths;"
  echo "[update_assembly_list] REMOVED species are left in the list untouched either way)."
  rm -f "$CUR_SORTED"
  exit 0
fi

# Apply: keep REMOVED species' rows as-is (don't delete data for species
# that merely vanished from this listing pass), update CHANGED paths, add
# NEW species.
awk -F'\t' 'NR==FNR{fresh[$1]=$2; next} {if ($1 in fresh) print $1"\t"fresh[$1]; else print}' \
  "$FRESH" "$ASSEMBLY_LIST" > "${ASSEMBLY_LIST}.tmp"
awk -F'\t' 'NR==FNR{seen[$1]=1; next} !($1 in seen)' "$ASSEMBLY_LIST" "$FRESH" >> "${ASSEMBLY_LIST}.tmp"
sort -o "${ASSEMBLY_LIST}.tmp" "${ASSEMBLY_LIST}.tmp"
mv "${ASSEMBLY_LIST}.tmp" "$ASSEMBLY_LIST"
echo
echo "[update_assembly_list] Applied: $ASSEMBLY_LIST now has $(wc -l < "$ASSEMBLY_LIST") species."

if [[ "$PURGE_STALE" == "1" && "$n_changed" -gt 0 ]]; then
  echo "[update_assembly_list] Purging stale per-species outputs for CHANGED species..."
  while read -r species; do
    [[ -z "$species" ]] && continue
    rm -fv \
      "${OUTPUT_DIR}/tr_tbls/${species}.tbl" \
      "${OUTPUT_DIR}/tr_tbls/${species}_all_tr_seqs_020226.tbl" \
      "${OUTPUT_DIR}/tert_tbls/${species}"*.tbl \
      "${OUTPUT_DIR}/tert_tblastn/${species}.tbl" \
      "${OUTPUT_DIR}/tidk/${species}.tidk.tsv" \
      "${OUTPUT_DIR}/tr_domains/${species}.tsv" \
      "${OUTPUT_DIR}/tr_template_loci/${species}.tsv" \
      2>/dev/null || true
  done <<< "$CHANGED"
  echo "[update_assembly_list] Re-run the per-stage scripts (run_tr_hmms.bash, run_tidk.bash,"
  echo "[update_assembly_list] run_tblastn_tert_bsub.bash, run_extract_tr_domains.bash, ...) to regenerate."
fi

rm -f "$CUR_SORTED"
