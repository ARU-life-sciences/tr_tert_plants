#!/usr/bin/env bash
set -euo pipefail

# Downloads reference genomes for species outside the DToL assembly list
# (inputs/dtol_plant_paths.txt) that are useful as extra tblastn targets --
# e.g. close relatives of a species with a suspicious TERT result, to sanity
# check whether a rice-only (or otherwise distant) query set is failing to
# find a real, divergent TERT rather than TERT being genuinely absent.
#
# Downloaded assemblies are NOT committed to git (see .gitignore); this
# script re-fetches them from NCBI. Each entry is recorded in
# inputs/extra_genomes.txt (species<TAB>path), which uses the same
# two-column format as dtol_plant_paths.txt, so it can be pointed at
# run_tblastn_tert_bsub.bash (ASSEMBLY_LIST=...) or run_tblastn_tert_local.bash.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
GENOME_DIR="${INPUT_DIR}/extra_genomes"
MANIFEST="${INPUT_DIR}/extra_genomes.txt"

mkdir -p "$GENOME_DIR"
touch "$MANIFEST"

# species, GCA/GCF accession, assembly name (as used in the NCBI FTP path)
entries=(
  "Santalum_album GCA_043873815.1 ASM4387381v1"
)

for entry in "${entries[@]}"; do
  read -r species accession asm_name <<< "$entry"

  species_dir="${GENOME_DIR}/${species}"
  out_gz="${species_dir}/${accession}.fna.gz"

  if [[ "${FORCE:-0}" != "1" ]] && [[ -s "$out_gz" ]]; then
    echo "[fetch_extra_genomes] ${species}: $out_gz exists; skipping. Set FORCE=1 to re-fetch."
  else
    mkdir -p "$species_dir"
    # NCBI FTP path is derived from the accession: GCA_043873815.1 -> GCA/043/873/815
    a="${accession%%.*}"
    p1="${a:4:3}"; p2="${a:7:3}"; p3="${a:10:3}"
    url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${a:0:3}/${p1}/${p2}/${p3}/${accession}_${asm_name}/${accession}_${asm_name}_genomic.fna.gz"

    echo "[fetch_extra_genomes] ${species}: downloading $url"
    curl -sS -o "$out_gz" "$url"
  fi

  if ! grep -qF -- "${species}"$'\t'"${out_gz}" "$MANIFEST" 2>/dev/null; then
    printf '%s\t%s\n' "$species" "$out_gz" >> "$MANIFEST"
  fi
done

echo "[fetch_extra_genomes] Done. Manifest: $MANIFEST"
