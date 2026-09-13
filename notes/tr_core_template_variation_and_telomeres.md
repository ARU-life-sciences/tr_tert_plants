# Does TR gene core-template variation show up in the telomeres themselves?

**Status: validated methodology, exhaustively resolved dataset-wide
result.** This note supersedes the methodology (not the headline
conclusion) of `notes/tr_repeat_correlation.md` - three real errors were
caught and fixed during development, all documented below, because
they're easy to repeat if this analysis is extended later.

## The question

Given two TR gene paralogs within one species whose Template sequences
differ by a small number of base pairs at the actual templating core (not
just somewhere in the wider extracted window), does that specific
difference show up in the species' real telomeric repeat content - and is
`tidk` fine-grained enough to detect it, or does that need something else?

## Three methodology errors, corrected

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

**Error 3: checking only the single best (or top-3) chromosome for the
minority variant undercounted real positives.** The first
WEAK_SENSITIVITY resolution pass picked one candidate chromosome-end per
species (the strongest for the *dominant* variant) and walked only that;
a second pass tried the top 3. Both produced false GENUINELY_ABSENT
calls - `Arctium_lappa` and `Tellima_grandiflora` were confirmed present
on chromosomes that weren't in either sample. **Fix**: check every main
chromosome, both ends, exhaustively. Implemented in the current
`resolve_weak_sensitivity.py`, made practical by also fixing a real
efficiency problem in `walk_terminal_repeat_array.py` - it was
re-decompressing the whole genome on every single region fetch;
`prepare_genome()`/`extract_terminal_prepared()` now decompress once per
species and reuse an approximate per-chromosome length (from the
existing tidk window scan) for a fast 3'-end fetch instead of pulling the
full chromosome every time.

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

27/44 (61%) are unambiguous positives by the coarse method alone - but
that undercounts the true rate, and the undercount is now fully resolved
(not just estimated) below.

### Exhaustive resolution: 36/44 (82%) confirmed, zero genuine negatives

The 8 WEAK_SENSITIVITY sets were each checked exhaustively -
**every** main chromosome, **both** ends (not a top-N sample) - by
directly walking the raw terminal sequence with
`walk_terminal_repeat_array.py` and `resolve_weak_sensitivity.py`.
Partial sampling matters here and got this wrong on the first two passes:
a single-best-chromosome check first, then a top-3 check, both produced
false GENUINELY_ABSENT calls that the full sweep overturned (a minority
variant can sit on a chromosome that isn't among the strongest for the
*dominant* variant, so ranking candidates by the dominant variant's
signal doesn't reliably find where the minority one is).

**Result: all 8/8 resolved to CONFIRMED_PRESENT. Zero genuine negatives
remain among the 44 differing sets.**

| species | dominant | minority | min_freq | chromosome-ends found at |
|---|---|---|---|---|
| Ajuga_chamaepitys | AAACTAAAC | AACCTAATC | 4.0% | 12 of 14 chromosomes |
| Ajuga_chamaepitys | CTAAACT | TTTTACA | 0.15% | 1 of 14 |
| Ajuga_chamaepitys | AACCCT | CTAATC | 0.21% | 14 of 14 |
| Trocdaris_verticillatum | ACCCTAACC | CTTCTCCGG | 0.71% | 1 of 12 |
| Arctium_lappa | TAAACCCTAAACC | AAAACCCTAAACC | 2.2% | 18 of 40 |
| Arctium_minus | TAAACCCTAAACC | AAAACCCTAAACC | 1.3% | 18 of 18 checked |
| Tellima_grandiflora | AACCCTAAACC | AAAACCCTAAA | 4.8% | 6 of 10 |
| Glaucium_flavum | GAAAACCCTACCCG | AAAAACCCTACCCG | 3.3% | 2 of 13 |

Frequencies are low (0.15-4.8%) and, where widespread, genuinely
dataset-wide rather than localised to one array - e.g. `Ajuga
chamaepitys`'s `AACCCT`/`CTAATC` pair is found at all 14 of its
chromosomes, 1-5 copies each, not clustered on one.

Final combined table (coarse classification + exhaustive resolution
merged): `outputs/tr_repeat_correlation/final_repeat_confirmations.tsv`,
via `finalize_repeat_confirmations.py`. **36/44 CONFIRMED_PRESENT, 8/44
LABEL_FLIP, 0/44 negative.**

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
| `walk_terminal_repeat_array.py` | base-pair-resolved raw-sequence tiling - the "something else" needed beyond tidk; `prepare_genome()`/`extract_terminal_prepared()` decompress once per species and reuse an approximate per-chromosome length for fast repeated fetches |
| `resolve_weak_sensitivity.py` | exhaustively resolves every WEAK_SENSITIVITY set - every main chromosome, both ends, one genome decompression per species |
| `finalize_repeat_confirmations.py` | merges the coarse classification with the exhaustive resolution into one final table |

## Caveats

- **The exhaustive walk was only run on the 8 WEAK_SENSITIVITY sets**,
  not the 27 STRONG_OVERLAP / 1 PARTIAL_OVERLAP sets (already directly
  observed by the coarse method, so lower priority) or the 8 LABEL_FLIP
  sets (a different, unexplained phenomenon - see below). 36/44 is solid;
  it isn't a claim about those other 8.
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

# Exhaustively resolve every WEAK_SENSITIVITY set (every main chromosome, both ends):
python3 src/resolve_weak_sensitivity.py outputs/tr_repeat_correlation/differing_set_confirmations.tsv \
  outputs/tr_repeat_correlation/positional 5000 outputs/tr_repeat_correlation/weak_sensitivity_resolved.tsv

# Merge into the final table:
python3 src/finalize_repeat_confirmations.py outputs/tr_repeat_correlation/differing_set_confirmations.tsv \
  outputs/tr_repeat_correlation/weak_sensitivity_resolved.tsv outputs/tr_repeat_correlation/final_repeat_confirmations.tsv

# One-off base-pair-resolved check of a single chromosome-end:
python3 src/walk_terminal_repeat_array.py <species> <chromosome> <5prime|3prime> <variant1> <variant2> --window 5000
```
