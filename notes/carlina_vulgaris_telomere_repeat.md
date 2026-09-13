# Carlina vulgaris: candidate non-canonical telomere repeat

**Status: candidate finding, not confirmed.** Computational evidence from two
independent methods converges on the same conclusion, but there's no
orthogonal (e.g. literature, wet-lab) confirmation yet.

## Summary

*Carlina vulgaris* (carline thistle, Asteraceae) appears to have replaced the
canonical land-plant telomere repeat (`TTTAGGG`) with a derived ~6 bp repeat,
`AAACTG`-family, at most (6 of 7) of its TR-homologous genomic loci. This
showed up first as an anomaly in automated TR-gene Template-domain
extraction (see `src/extract_tr_domains.py`), then was independently
corroborated by TIDK positional search directly on the genome.

Of the 420 DToL species with a confident TR gene hit, this is the only one
(after correcting a scoring bug that had also flagged a couple of others -
see below) where the extracted Template region shows no trace of the
canonical repeat.

## Background

The TR gene (telomerase RNA) contains a short template region that TERT
reverse-transcribes into the telomeric repeat. `src/extract_tr_domains.py`
locates this template computationally for every DToL species (see
`inputs/TR_domains.all_tr_seqs.tsv` for how the domain is anchored) and
cross-checks the extracted sequence's reverse complement against both the
canonical repeat and that species' own TIDK-discovered repeat
(`outputs/tidk/*.tidk.tsv`). Almost all species (415/420, 98.8%) confirm
canonical `TTTAGGG` cleanly. Five did not on the first pass; two of those
turned out to be a matching-algorithm bug (picking the longest coincidental
match across candidates instead of prioritising the canonical repeat - fixed
in `extract_tr_domains.py`), and a third (*Stachys alpina*) is very likely
just an extraction-window registration issue, since its genome-wide TIDK
data shows an abundant canonical-family repeat as the #1 hit. That leaves
*Carlina vulgaris* (and a weaker, unresolved case in *Jasione montana* that
turned out on closer inspection to be ordinary point-level sequence
divergence around an otherwise correctly-registered canonical template, not
a different repeat family).

## Evidence

### 1. TR gene: divergent Template at most (not quite all) genomic loci

**Updated** after running the systematic within-species multi-locus survey
(`src/extract_tr_template_all_loci.py`, see
`notes/within_species_tr_template_variation.md`) - Carlina actually has
**7** TR-homologous loci passing the same E<=1e-4 threshold used dataset-wide
(not 4; the initial manual check found the strongest 4 and missed 3 weaker
ones), and **one of them is not part of the AAACTG pattern**:

```
outputs/tr_template_loci/Carlina_vulgaris.tsv
SUPER_2 (best, E=3.3e-19)  ACAAACTGAAACTGCCTCTTTGGGTTCGATCCTTGGGTTC
SUPER_5 (E=7.2e-18)        TTAAAACCGAACTGTCCTCTGTGGGTTGCTTCTTGGGTTTC
SUPER_9 (E=2.0e-16)        TAACACCTTAAACCGAACTGTCCAATGGGTTGCTTCTTGGGTTTC
SUPER_5 (E=1.9e-15)        TTAAACCGAACTGTCCTCTGTGGGTTGCTTGTTGGGTTTC
SUPER_2 (E=1.8e-11)        CAAAACTAAACTGGTCTCTTATGAGGGTGTGGTGCTTGGGTTTC
SUPER_1 (E=3.8e-8)         AAACTGTCCTCTTTGGGTTCGATCCTTGGGTTC
SUPER_4 (E=9.7e-6)         ACAAACCCTAACCTTTCTTAGGGTGTGTCTTGGGGTCA   <- canonical-matching
```

Six loci (E = 2e-16 to 3.8e-8, i.e. the strongest, most confident hits) carry
the `GAACTG`/`AAACTG` motif and don't resemble canonical (compare to *Dryas
octopetala*, a confirmed canonical control, whose Template reads
`ATTAACCCTAAACCCTCCCATT-G-CGGGATTT---ATGTGGGGGTC`). The seventh, **weakest**
locus (SUPER_4, E=9.7e-6 - roughly four orders of magnitude less significant
than the rest) reverse-complements cleanly to canonical `TTTAGGG`
(`AAACCCTAA` found via `find_repeat_in_window`). So this isn't a clean
"every copy switched" story - it's the same **mixed** pattern seen at 126 of
331 (38%) multi-locus species genome-wide (see the within-species note):
most loci diverged, one weak one didn't. Given SUPER_4 is also the least
significant hit of the seven, it's plausible it's either a genuinely
ancestral/relict copy that hasn't been through whatever process altered the
others, or a marginal HMM hit that isn't reliably homologous at all - this
note doesn't resolve which.

None of this changes the strongest evidence (section 2, below): the flanking
domains (G-rich 5' end, Conserved region) look normal and homologous across
all loci, so these are genuinely TR-gene-like hits, not off-target matches,
and the `AAACTG` motif shared by six of seven independent loci also happens
to match Carlina's own #2 most abundant genome-wide repeat by TIDK's
`explore` mode (`outputs/tidk/Carlina_vulgaris.tidk.tsv`:
`AAACTGAAACTG...`, 12,181 runs).

### 2. Genome: real chromosome-terminal enrichment for AAACTG, none for TTTAGGG

`tidk explore`'s output ranks candidates by count within a proportion of
chromosome length from each end (default `--distance 0.01`), so its "most
abundant" ranking is already terminal-proximal, not literally genome-wide -
but to check the actual positional pattern (rather than relying on that
distance-restricted count as if it were an unrestricted spatial argument),
`tidk search` was run directly with 5 kb windows across the primary
assembly (`daCarVulg1.hap1.2.primary.fa.gz`) for both `AAACTG` and
`TTTAGGG`:

```bash
tidk search -s AAACTG  -w 5000 -o carlina_AAACTG  -d <dir> --log <fasta>
tidk search -s TTTAGGG -w 5000 -o carlina_TTTAGGG -d <dir> --log <fasta>
```

`AAACTG` shows a textbook telomere signature - forward-strand repeat counts
spike at the very first window(s) of several chromosomes, reverse-strand
counts spike at the very last window(s) of others, decaying to ~0 within
the first few windows moving inward:

| chromosome | terminal window counts | end |
|---|---|---|
| SUPER_2 | 90-171 | 5' |
| SUPER_5 | 183-237 | 5' |
| SUPER_7 | 191-322 | 5' |
| SUPER_3 | 105-179 (reverse strand) | 3' |
| SUPER_4 | 118-162 (reverse strand) | 3' |
| SUPER_8 | 114-186 (reverse strand) | 3' |
| SUPER_9 | 85-187 (reverse strand) | 3' |

`TTTAGGG` shows none of this at any of the 10 main chromosomes (SUPER_1-10)
- counts are ~0 at every terminus, with only weak, scattered interior peaks
(max ~20-32, not chromosome-proximal). This is a meaningful absence, not
just lower abundance.

## Interpretation

Two independent lines of evidence (TR gene homology + genome-wide
positional repeat search) converge: this species' actual telomeric repeat
looks like an `AAACTG`-family unit, and most of its TR-homologous loci are
consistent with encoding that repeat rather than the canonical one (one weak
locus is the exception - see above). This would be the same class of event
documented in Fajkus et al. 2019 (NAR,
[10.1093/nar/gkz695](https://doi.org/10.1093/nar/gkz695)) for Allium,
Cestrum and Scilla - independent, lineage-specific replacement of the
ancestral land-plant telomere repeat - just not a previously known instance,
and (based on the other 10 Asteraceae in this dataset, all canonical)
apparently restricted to this one species rather than shared family-wide.

## Caveats

- The Template-domain match to `AAACTG` isn't a clean textbook doubled-repeat
  match the way canonical hits are elsewhere in this dataset - it shows up
  as a single copy (occasionally a tandem pair at the best-scoring locus)
  rather than the ~1.5 repeat-unit span you'd expect from a precisely
  registered template. The automated matcher in `extract_tr_domains.py` is
  tuned for the canonical 7 bp case and may not be finding the true
  boundaries of a shorter, possibly degenerate 6 bp unit precisely.
- All seven TR-homologous loci could be a paralogous/duplicated gene family
  rather than independent confirmation of a single functional gene - i.e.
  these are correlated data points from a shared ancestry, not seven fully
  independent lines of evidence, though the six-vs-one split does rule out
  single-locus assembly or extraction artifacts as the explanation.
- Which locus (if any) is the actual functionally active TR gene is unknown -
  it could be one of the six AAACTG-pattern loci, the one canonical-looking
  locus, or none of them (i.e. the real functional copy might be elsewhere
  and simply not detected).
- No literature check has been done to see whether this is already a known
  finding for this species or genus.
- Assembly: `daCarVulg1.hap1.2.primary.fa.gz` (haplotype 1, curation round
  2) - worth confirming the second haplotype agrees, if available.

## Reproducing

```bash
# Domain map + per-species (best-hit only) extraction:
python3 src/map_tr_domains.py outputs/tr_tbls/all_tr_seqs.hmm inputs/TR_domains.all_tr_seqs.tsv
bash src/run_extract_tr_domains.bash
grep Carlina_vulgaris outputs/tr_domains/all_species_tr_domains.tsv

# All-loci extraction (what section 1 above is based on):
python3 src/extract_tr_template_all_loci.py Carlina_vulgaris \
  "$(awk -F'\t' '$1=="Carlina_vulgaris"{print $2}' inputs/dtol_plant_paths.txt)" \
  outputs/tr_tbls/Carlina_vulgaris_all_tr_seqs_020226.tbl \
  outputs/tr_tbls/all_tr_seqs.hmm \
  outputs/tr_template_loci

# TIDK positional search (not part of the standard pipeline - run ad hoc):
tidk search -s AAACTG  -w 5000 -o carlina_AAACTG  -d <outdir> --log \
  /data/tol/data/darwin/dicots/Carlina_vulgaris/assembly/curated/daCarVulg1.2/daCarVulg1.hap1.2.primary.fa.gz
tidk search -s TTTAGGG -w 5000 -o carlina_TTTAGGG -d <outdir> --log \
  /data/tol/data/darwin/dicots/Carlina_vulgaris/assembly/curated/daCarVulg1.2/daCarVulg1.hap1.2.primary.fa.gz
```
