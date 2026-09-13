# TR/TERT Plant Genome Scanner

A bioinformatics pipeline for identifying telomerase RNA (TR) genes, Telomerase Reverse Transcriptase (TERT) genes, and telomeric repeat motifs across plant genome assemblies from the Darwin Tree of Life (DToL) project.

## Overview

The pipeline processes ~488 plant species assemblies and uses three complementary approaches:

- **TR detection** — HMM profile search (nhmmer)
- **TERT detection** — HMM profile search + protein BLAST (tblastn)
- **Telomeric repeats** — tandem repeat discovery (TIDK)

## Dependencies

- [HMMER](http://hmmer.org/) 3.4 (`nhmmer`)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) 2.16.0+ (`tblastn`)
- [TIDK](https://github.com/tolkit/telomeric-identifier)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (for rebuilding HMM profiles)
- [CD-HIT](https://sites.google.com/view/cd-hit) (for TERT input clustering)
- LSF batch system (`bsub`)

## Repository Structure

```
inputs/        # HMM profiles, query sequences, assembly list
outputs/
  tr_tbls/     # TR HMM search results
  tr_domains/  # Per-species TR gene domain extraction (USE, TATA-box, Template, ...)
  tert_tbls/   # TERT HMM search results
  tert_tblastn/ # TERT BLAST results
  tidk/        # Telomeric repeat results
src/           # Pipeline scripts
notes/         # Write-ups of individual findings that came out of the pipeline
```

## Usage

Scripts in `src/` are run in approximate order:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `fetch_latest_assemblies.bash` | Retrieve DToL plant assembly paths |
| 2 | `fetch_tert.bash` | Download reference TERT sequences from NCBI |
| 3 | `run_tr_hmms.bash` | Search assemblies for TR genes |
| 4 | `get_all_tr_seqs.bash` | Extract and consolidate TR hits |
| 5 | `recreate_TR_hmm.bash` | Rebuild TR HMM from identified sequences |
| 6 | `run_tert_hmms.bash` | Search assemblies for TERT (3 HMM profiles) |
| 7 | `run_tblastn_tert_bsub.bash` | Submit LSF jobs for TERT protein BLAST |
| 8 | `run_tidk.bash` | Discover telomeric repeat motifs |
| 9 | `map_tr_domains.py` | Map TR gene domains (USE, TATA-box, Template, ...) onto HMM columns |
| 10 | `run_extract_tr_domains.bash` | Extract each domain per species, cross-check Template against TIDK |
| 11 | `run_extract_tr_template_all_loci.bash` | Extract Template from every TR-homologous locus per species (within-species variation) |
| 12 | `tr_repeat_correlation.py` | Screen every multi-locus species for Template-vs-telomere-repeat sequence correlation |
| 13 | `tr_repeat_correlation_triage.py` | Shortlist the strongest candidates from step 12 for positional confirmation |
| 14 | `confirm_repeat_positional.py` | Confirm a candidate repeat is genuinely chromosome-terminal via `tidk search` |
| 15 | `tr_core_template_diff.py` | Position-anchored pairwise diff of TR loci's true templating core (not the full padded window) |
| 16 | `tr_core_template_sets.py` | Species-level (not pairwise) summary of core-template identity/variation across TR gene copies |
| 17 | `run_confirm_differing_sets.bash` / `_bsub.bash` | Positionally confirm every differing set's dominant + minority variant |
| 18 | `summarize_repeat_confirmations.py` | Classify each differing set's telomeric confirmation result |
| 19 | `walk_terminal_repeat_array.py` | Base-pair-resolved raw-sequence check for variants too rare for `tidk search` to see |

## Onboarding new genomes

DToL adds and re-curates plant assemblies over time, so `inputs/dtol_plant_paths.txt`
goes stale. To pick up new/updated genomes:

```bash
# 1. See what's changed (report-only, doesn't touch anything):
bash src/update_assembly_list.bash

# 2. Apply it once you're happy with the diff:
bash src/update_assembly_list.bash --apply
#    add --purge-stale too if you want CHANGED species' existing outputs
#    (tr_tbls, tert_tblastn, tidk, tr_domains, tr_template_loci) deleted so
#    they get regenerated against the new assembly, rather than silently
#    keeping stale coordinates from the old one.

# 3. Check what's left to run:
bash src/pipeline_status.bash

# 4. Run whichever per-stage scripts are needed (steps 3-11 above) - each
#    one already skips species with existing output, so it's safe to just
#    re-run them after step 1/2 and they'll only process what's new.
```

`update_assembly_list.bash` reports three categories: **NEW** species (not
seen before - the common case), **CHANGED** (the recorded "latest" assembly
path differs from what's on disk now - this happens on re-curation, and can
mean genuinely different sequence content under an unchanged filename, not
just a rename: see `notes/`), and **REMOVED** (dropped out of the assembly
tree entirely - never acted on automatically).

Note that `run_tblastn_tert_bsub.bash` in particular can be very expensive
per genome (one real case took 12+ days for a single large assembly) - for
a big batch of new species, expect that stage to dominate wall-clock time by
a wide margin over everything else in the pipeline.

## Findings

Notable results from the pipeline are written up in `notes/`:

- [`carlina_vulgaris_telomere_repeat.md`](notes/carlina_vulgaris_telomere_repeat.md) - candidate non-canonical (`AAACTG`-family) telomere repeat
- [`within_species_tr_template_variation.md`](notes/within_species_tr_template_variation.md) - 38% of multi-locus species carry both canonical- and divergent-Template loci in the same genome
- [`organelle_and_cobiont_assembly_pitfall.md`](notes/organelle_and_cobiont_assembly_pitfall.md) - two species' "latest assembly" silently pointed at a plastid/cobiont genome instead of the nuclear one; now screened automatically
- [`tr_repeat_correlation.md`](notes/tr_repeat_correlation.md) - fine-scale TR-paralog point mutations track specific chromosome-terminal repeat variants, confirmed in 8/8 species checked (generalizes the Carlina finding dataset-wide)
- [`tr_core_template_variation_and_telomeres.md`](notes/tr_core_template_variation_and_telomeres.md) - properly position-anchored core-template comparison: 71-84% of TR gene copies share an identical core even at high copy number, and most (>=27/44, likely more) of the differing cases show the specific point-mutant variant genuinely present at the telomere - includes two real methodology corrections worth reading before extending this analysis
