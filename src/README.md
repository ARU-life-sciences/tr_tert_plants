## Pipeline overview

This project scans plant genome assemblies for:

- **TR (telomerase RNA) genes**
- **TERT (telomerase reverse transcriptase) genes**
- **Telomeric repeat motifs**

It assumes this layout:

- `inputs/`  – query sequences, HMMs, assembly list
- `outputs/` – all results
- `src/`     – pipeline scripts

Assemblies are listed in:

- `inputs/dtol_plant_paths.txt` (format: `species<TAB>/path/to/assembly.fa[.gz]`)

---

## Scripts (what they do)

### 0. Assemblies

**`fetch_latest_assemblies.bash`**  
Gets latest DToL assembly paths.

- `inputs/dtol_plant_paths.txt`

---

### 1. TERT query sequences

**`fetch_tert.bash`**  
Fetches reference TERT sequences:

- `inputs/plant_TERT_regions.fasta`  – nucleotide regions for TERT
- `inputs/plant_TERT_regions.faa`    – protein sequences for TERT

Skips if both outputs already exist (unless `FORCE=1`).

---

### 2. TR HMM search

**`run_tr_hmms.bash`**  
Runs `nhmmer` with `inputs/TR.hmm` against each assembly.

Output per species:

- `outputs/tr_tbls/<species>.tbl`

Skips species with existing `.tbl`.

---

### 3. Extract TR sequences

**`get_all_tr_seqs.bash`**  
Parses all `outputs/tr_tbls/*.tbl` with `extract_nhmmer_tblout` and concatenates hits:

- `outputs/tr_tbls/all_tr_seqs.fasta`

Used to realign TRs or rebuild the TR HMM.

---

### 4. Rebuild TR HMM (optional)

**`recreate_TR_hmm.bash`**  
From `all_tr_seqs.fasta`, runs MAFFT + HMMER:

- `outputs/tr_tbls/all_tr_seqs.aln`
- `outputs/tr_tbls/all_tr_seqs.hmm` (+ pressed `.h3*` files)

Skips if HMM already exists (unless `FORCE=1`).

---

### 5. TERT HMM search (nhmmer)

**`run_tert_hmms.bash`**  
Uses three TERT HMMs:

- `inputs/clusters/isoform_X1_tert.hmm`
- `inputs/clusters/catalytic_subunit_tert.hmm`
- `inputs/clusters/tert1.hmm`

Submits one LSF job per species, running all three `nhmmer` searches in series.

Outputs per species:

- `outputs/tert_tbls/<species>/<species>_isoformx1.tbl`
- `outputs/tert_tbls/<species>/<species}_catalytic_subunit.tbl`
- `outputs/tert_tbls/<species>/<species>_tert1.tbl`

Skips species where all three files exist and are non-empty.

---

### 6. TERT tblastn search

**`run_tblastn_tert_bsub.bash`**  
Submits LSF jobs for tblastn using TERT protein queries.

**`worker_tblastn_tert.bash`**  
Per species:

- builds a temporary BLAST DB from the assembly
- runs `tblastn` with `inputs/plant_TERT_regions.faa`

Output per species:

- `outputs/tert_tblastn/<species>.tbl`

Skips if `.tbl` already exists.

---

### 6b. Extra (non-DToL) comparison genomes

**`fetch_extra_genomes.bash`**
Downloads reference genomes for species that aren't in the DToL assembly list but are
useful tblastn targets — e.g. a close relative to sanity-check a suspicious TERT
result against a genuinely divergent query. Writes to `inputs/extra_genomes/<species>/`
(gitignored — re-run the script rather than committing genomes) and records each one in:

- `inputs/extra_genomes.txt` (same two-column format as `dtol_plant_paths.txt`)

**`run_tblastn_tert_local.bash`**
Non-LSF counterpart to `run_tblastn_tert_bsub.bash`: runs `worker_tblastn_tert.bash`
directly, in series, for every line in `inputs/extra_genomes.txt` (or any assembly
list passed via `ASSEMBLY_LIST`). Falls back to `tblastn`/`makeblastdb` on `$PATH`
if the Sanger HPC `BLASTBIN` default isn't present. Output lands in the same place
as the bsub path:

- `outputs/tert_tblastn/<species>.tbl`

Currently populated with **`Santalum_album`** (GCA_043873815.1, T2T reference) —
added because Santalales has no annotated TERT/telomerase records anywhere in
GenBank (checked order-wide via NCBI eutils), so the rice-only query set in
`plant_TERT_regions.faa` is a weak negative test for TERT absence in *Viscum*.
`plant_TERT_regions.faa` also gained a manually-curated query,
`TERT_Arceuthobium_sichuanense_exon9` (a partial exon-9 fragment recovered by
BLAST from raw reads, BioProject PRJNA307530 — not from GenBank, see
`fetch_tert.bash`), to bracket *Viscum* with a closer relative on both sides
(Viscaceae via Arceuthobium, wider Santalaceae via Santalum).

**`extract_hit_region.py`**
Given a genome, a chrom:start-end/strand, and a name, extracts the region and
translates it in whichever frame has the fewest stop codons (stderr reports
all three, so the choice is auditable) — for turning an unannotated tblastn
hit into a new curated query. Used to pull the *Santalum album* TERT locus
found by `run_tblastn_tert_local.bash` (chr1 `CP159901.1:20959560-20960756`,
minus strand — the single strongest, near gap-free hit block across the rice
query set, e-value down to 9.4e-110): translates clean with zero stop codons
in frame 0, and shares a clearly homologous motif with the independently
obtained `TERT_Arceuthobium_sichuanense_exon9`
(`...ANF[FY]VTES[EY]G[QK]K[IT]VYYYRKS[VI]WE[KN]L...`), which cross-validates
both as genuine TERT despite neither being annotated anywhere. Added to
`plant_TERT_regions.fasta`/`.faa` as `TERT_Santalum_album_chr1locus`:

```
python3 src/extract_hit_region.py \
  inputs/extra_genomes/Santalum_album/GCA_043873815.1.fna.gz \
  CP159901.1 20959560 20960756 - TERT_Santalum_album_chr1locus
```

---

### 7. Telomeric repeat discovery

**`run_tidk.bash`**  
Runs `tidk explore` on each assembly to find abundant tandem repeats near chromosome ends.

Output per species:

- `outputs/tidk/<species>.tidk.tsv`

Contains candidate telomeric repeat units and counts.  
Skips if `.tidk.tsv` exists (unless `FORCE=1`).

---

## Recommended run order

From the project root:

1. `bash src/fetch_latest_assemblies.bash`
2. `bash src/fetch_tert.bash`
3. `bash src/run_tr_hmms.bash`
4. `bash src/get_all_tr_seqs.bash`
5. *(optional)* `bash src/recreate_TR_hmm.bash`
6. `bash src/run_tert_hmms.bash`
7. `bash src/run_tblastn_tert_bsub.bash`
7b. *(optional)* `bash src/fetch_extra_genomes.bash && bash src/run_tblastn_tert_local.bash`
8. `bash src/run_tidk.bash`

After this, you have for each species:

- TR hits: `outputs/tr_tbls/<species>.tbl`
- Combined TR sequences: `outputs/tr_tbls/all_tr_seqs.fasta`
- TERT HMM hits: `outputs/tert_tbls/<species>/*.tbl`
- TERT tblastn hits: `outputs/tert_tblastn/<species>.tbl`
- Telomeric repeats: `outputs/tidk/<species>.tidk.tsv`

