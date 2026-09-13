# Does TR gene core-template variation show up in the telomeres themselves?

**Status: validated methodology, systematic dataset-wide result.** This
note supersedes the methodology (not the headline conclusion) of
`notes/tr_repeat_correlation.md` - two real errors were caught and fixed
during development, both documented below, because they're easy to
repeat if this analysis is extended later.

## The question

Given two TR gene paralogs within one species whose Template sequences
differ by a small number of base pairs at the actual templating core (not
just somewhere in the wider extracted window), does that specific
difference show up in the species' real telomeric repeat content - and is
`tidk` fine-grained enough to detect it, or does that need something else?

## Two methodology errors, corrected

Both were caught by direct questioning of an initial result, not by
internal review - worth reading before extending this analysis.

**Error 1: "core template" was initially the whole padded extraction
window, not the true templating region.** `outputs/tr_template_loci`
extracts a fixed 47-HMM-column window (chosen wide, to safely cover
template length variation across ~330 species) - most of that window is
flanking scaffold sequence around a much shorter real repeat-encoding
core. An initial full-window Hamming-distance comparison between two
*Cornus sanguinea* loci found 8 differences and wrongly implied the
templating core itself differed substantially; checking exactly where
those 8 differences fell showed 7 of them were outside the actual
matched-repeat region entirely. **Fix**: define "the actual biological
template" as the locus's own best matched chunk (the region
independently verified to be shared with a real genome repeat, from
`tr_repeat_correlation.py`), anchored to its offset within the extracted
window - and only compare two loci's cores when they anchor at the same
(+/-1) offset, i.e. genuinely the same structural position, not a
same-length coincidence from unrelated parts of the window. Implemented
in `tr_core_template_diff.py`.

**Error 2: the "97% of same-species locus pairs share an identical core"
statistic was inflated by pairwise counting.** A species with N mutually
identical loci contributes C(N,2) trivially-identical pairs - a single
44-copy fully-conserved species (`Levisticum_officinale`) alone
contributes 946 such pairs, swamping genuinely variable species in the
aggregate. **Fix**: group loci per species into comparable sets (same
core length, same anchor position) and report set-level, not pair-level,
statistics. Implemented in `tr_core_template_sets.py`.

## Corrected findings

### How often do TR gene copies actually differ at the true core?

267 comparable (species, core-length) sets with >=2 loci:

| | sets | % |
|---|---|---|
| All loci share an identical core | 223 | 83.5% |
| Restricted to sets with >=5 loci (38 sets) | 27 | 71.1% |

**Higher copy number does not mean more conservation - if anything the
opposite.** More copies means more chances for a difference to arise and
persist, and the data bears that out: high-copy sets are *less* uniformly
identical (71.1%) than the full set (83.5%). Concrete examples:

| species | core length | loci | distinct cores |
|---|---|---|---|
| Levisticum_officinale | 7bp | 44 | 1 (perfectly conserved) |
| Anagallis_arvensis | 6bp | 17 | 6 |
| Solanum_nigrum | 11bp | 12 | 3 |
| Empetrum_nigrum | 12bp | 12 | 1 |
| Lythrum_salicaria | 11bp | 6 | 3 (half the copies differ) |

### When cores do differ, does the telomere reflect it?

44 differing sets (>1 distinct core) across 33 species. For each, ran
`confirm_repeat_positional.py` (window-count `tidk search` threshold) on
both the dominant (most TR-copies) and top minority variant, classified
via `summarize_repeat_confirmations.py`:

| classification | sets | meaning |
|---|---|---|
| STRONG_OVERLAP | 27 | both variants enriched, at largely the same chromosomes |
| PARTIAL_OVERLAP | 1 | both enriched, less shared overlap |
| LABEL_FLIP | 8 | the by-TR-copy-count "dominant" variant shows **no** telomeric signal; the "minority" one is the genome's real dominant repeat |
| WEAK_SENSITIVITY | 8 | dominant enriched, minority shows no/weak signal by this coarse method |

**27/44 (61%) are unambiguous positives** by the coarse method alone. But
that undercounts the true rate - see the calibration result below.

### Calibration: the coarse method misses real, low-frequency variants

Two WEAK_SENSITIVITY-by-coarse-method cases (`Cornus_sanguinea`,
`Arctium_lappa`) were independently checked by directly walking the raw
terminal genomic sequence base-by-base
(`walk_terminal_repeat_array.py`, which tiles the actual sequence rather
than counting window occurrences) - **both showed the minority variant
genuinely present, interspersed throughout the array, at ~1-3% and ~12%
frequency respectively.** A real, low-frequency, interspersed variant
easily fails `confirm_repeat_positional.py`'s fold-over-median threshold
(designed to detect a variant abundant enough to dominate a window) while
still being unambiguously real. **So the WEAK_SENSITIVITY bucket is a mix
of genuine negatives and under-detected true positives, not confirmed
negatives** - distinguishing them requires the slower walk tool, which
wasn't run at full scale (see Caveats).

One clean genuine negative was also found this way for comparison:
*Ajuga chamaepitys*'s `AACCTAATC` variant (dominant `AAACTAAAC` enriched
at 10/14 chromosomes) shows zero signal even by direct positional search
- a real absence, not a sensitivity artifact, since it's a completely
unrelated sequence rather than a 1bp neighbour of something already
confirmed present.

### Unexpected side finding: TR-copy count doesn't predict telomeric dominance

8/44 sets are LABEL_FLIP - the variant with *more* independent TR gene
copies backing it is telomerically silent, while the "minority" variant
(fewer TR copies) is the genome's actual dominant repeat (e.g.
*Tilia_cordata*: 0/40 vs 40/40 chromosomes). TR gene paralog count is not
a reliable proxy for which sequence a species' telomeres actually use.

## Scripts

| script | role |
|---|---|
| `tr_repeat_correlation.py` | (existing) per-locus best match to a genome repeat + shuffle-control significance |
| `tr_core_template_diff.py` | position-anchored pairwise core-template comparison (fixes Error 1) |
| `tr_core_template_sets.py` | species-level comparable-set grouping and identical/differing summary (fixes Error 2); also identifies dominant/minority variants per differing set |
| `run_confirm_differing_sets.bash` / `_bsub.bash` | run `confirm_repeat_positional.py` for every differing set's dominant + minority variant (local GNU-parallel or LSF) |
| `summarize_repeat_confirmations.py` | classify each differing set's confirmation result (STRONG_OVERLAP / PARTIAL_OVERLAP / LABEL_FLIP / WEAK_SENSITIVITY / NEITHER) |
| `walk_terminal_repeat_array.py` | base-pair-resolved raw-sequence tiling - the "something else" needed beyond tidk for the calibration check above; not run at full 44-set scale (see Caveats) |

## Caveats

- **The walk-tool calibration was only done for 2 of the 8
  WEAK_SENSITIVITY sets.** The true STRONG_OVERLAP rate is very likely
  higher than 27/44 (61%) once the rest are checked the same way, but
  that hasn't been done - don't quote 61% as a ceiling.
- **Haplotype-phasing risk is real, not hypothetical.** While building
  this, `Empetrum_nigrum`'s TR loci turned out to be spread across
  `SUPER_N_HAP2/HAP3/HAP4`-style scaffold names, consistent with a
  multi-haplotype-phased (possibly polyploid) assembly - meaning some
  "independent loci" could be the same physical locus's different
  haplotype copies (allelic variation) rather than independent
  paralogous gene copies. This wasn't systematically checked or
  filtered across the dataset; the `_HAP\d+` suffix pattern is a cheap
  thing to screen for before trusting a "many independent copies" claim
  for a specific species.
- **Orientation ambiguity remains unresolved** (see
  `tr_repeat_correlation.md`) - most core-template matches are in the
  Template's `direct` sense rather than the textbook `revcomp`
  templating-mechanism expectation, and this analysis didn't revisit
  that question.
- **LABEL_FLIP cases need their own explanation**, not yet investigated -
  is the "dominant by TR-copy-count" variant a decayed/pseudogenized
  paralog cluster unrelated to the actual functional gene, or something
  else? Open question.

## Reproducing

```bash
# Already have this from tr_repeat_correlation.md's stage 1:
#   outputs/tr_repeat_correlation/all_species.tsv

# Species-level comparable sets + differing-set identification:
python3 src/tr_core_template_sets.py outputs/tr_repeat_correlation/all_species.tsv \
  outputs/tr_repeat_correlation/core_template_sets.tsv

# Position-anchored pairwise view (used to find clean 1bp examples):
python3 src/tr_core_template_diff.py outputs/tr_repeat_correlation/all_species.tsv 1 \
  outputs/tr_repeat_correlation/core_template_diffs.tsv

# Confirm every differing set's dominant + minority variant positionally:
bash src/run_confirm_differing_sets.bash          # local
bash src/run_confirm_differing_sets_bsub.bash     # or via LSF

# Classify the results:
python3 src/summarize_repeat_confirmations.py outputs/tr_repeat_correlation/core_template_sets.tsv \
  outputs/tr_repeat_correlation/positional outputs/tr_repeat_correlation/differing_set_confirmations.tsv

# Base-pair-resolved calibration check for any WEAK_SENSITIVITY case:
python3 src/walk_terminal_repeat_array.py <species> <chromosome> <5prime|3prime> <variant1> <variant2> --window 5000
```
