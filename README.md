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
  tert_tbls/   # TERT HMM search results
  tert_tblastn/ # TERT BLAST results
  tidk/        # Telomeric repeat results
src/           # Pipeline scripts
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
