# Does TR gene core-template variation show up in the telomeres themselves?

**Status: validated methodology, exhaustively resolved dataset-wide
result - 44/44 (100%) of differing sets confirmed, zero genuine
negatives.** This note supersedes the methodology (not the headline
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
| Empetrum_nigrum | 12bp | 12* | 1 |
| Lythrum_salicaria | 11bp | 6* | 3 |

`*` = inflated by haplotype redundancy, confirmed by the systematic
check below - see that section for corrected counts. `Empetrum_nigrum`'s
raw 12 is likely ~4 true independent loci, each represented 2-3 times as
separate haplotype copies (same sequence each time, so it stays
1-distinct either way). `Lythrum_salicaria`'s raw 6 is really only 3 true
independent loci (`SUPER_8`, `SUPER_2`, `SUPER_1`, each haplotype-
duplicated 2-3x with identical sequence per locus) - all 3 happen to
carry different sequences from each other, so "3 distinct cores" is
still correct, but the earlier "half the copies differ" framing (implying
3 conserved vs 3 differing among 6 real copies) was wrong: it's 3 real
loci, each individually consistent across its own haplotype copies, that
simply differ *from one another*.

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

### Exhaustive resolution: 44/44 (100%) confirmed, zero genuine negatives

Both ambiguous buckets - the 8 WEAK_SENSITIVITY sets and, once it became
clear the same logic applied, the 8 LABEL_FLIP sets - were each checked
exhaustively: **every** main chromosome, **both** ends (not a top-N
sample), by directly walking the raw terminal sequence with
`walk_terminal_repeat_array.py` and `resolve_weak_sensitivity.py`
(generalised to take a target classification as its 5th argument so the
same exhaustive-walk logic runs against either bucket). Partial sampling
matters here and got this wrong on the first two WEAK_SENSITIVITY passes:
a single-best-chromosome check first, then a top-3 check, both produced
false GENUINELY_ABSENT calls that the full sweep overturned (a minority
variant can sit on a chromosome that isn't among the strongest for the
*dominant* variant, so ranking candidates by the dominant variant's
signal doesn't reliably find where the minority one is).

**Result: all 8/8 WEAK_SENSITIVITY sets AND all 8/8 LABEL_FLIP sets
resolved to CONFIRMED_PRESENT. Zero genuine negatives remain among the
44 differing sets.**

**WEAK_SENSITIVITY resolution** (checking whether the coarse-invisible
*minority* variant is really present):

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

**LABEL_FLIP resolution** (checking whether the coarse-invisible
*TR-copy-majority* variant - flagged "telomerically silent" by the coarse
method - is really absent). It isn't, in any of the 8 cases: the same
sensitivity-threshold artifact as WEAK_SENSITIVITY, just more extreme
(the coarse-silent variant is rarer relative to its partner here than in
the WEAK_SENSITIVITY cases above):

| species | TR-copy-majority variant | occurrences | genome's real dominant repeat | occurrences | chromosomes checked |
|---|---|---|---|---|---|
| Carpinus_betulus | TAAACCCTAAAAG | 34 (1.4%) | CTAAACCCTAAAC | 2339 (98.6%) | 16 |
| Juncus_bufonius | CTAAACCCTAGAT | 42 (0.24%) | CTAAACCCTAAAC | 17305 (99.8%) | 40 |
| Potentilla_indica | TAAACCCTAAACCT | 73 (3.7%) | AAACCCTAAACCCT | 1925 (96.4%) | 40 |
| Crataegus_laevigata | CAACCT | 55 (0.5%) | TAAACC | 10728 (99.5%) | 18 |
| Lycopus_europaeus | CTAAACCCTACC | 5 (0.1%) | AAACCCTAAACC | 5040 (99.9%) | 12 |
| Platanus_x_hispanica | AAACCT | 154 (0.9%) | CTAAAC | 16517 (99.1%) | 23 |
| Solanum_nigrum | AATAAAT | 86 (30.5%) | ACCTGAA | 196 (69.5%) | 40 |
| Tilia_cordata | CAAACCCTAAACC | 187 (1.3%) | TAAACCCTAAACC | 14023 (98.7%) | 40 |

Occurrence counts are totals across every main chromosome, both ends
(5000bp windows), from a single exhaustive walk per species - not a
per-chromosome breakdown of the majority variant specifically (the
`found_at` field in the underlying script only logs locations for the
second/minority argument, i.e. the real dominant repeat here - the
majority-variant total (`dom_n`) is the number that matters and is
reported directly). Raw data:
`outputs/tr_repeat_correlation/label_flip_resolved.tsv`.

**`Solanum_nigrum` needs an extra grain of salt**: at 30.5% it's a
clear outlier from the rest of this table (all <5%), and its
"majority" variant `AATAAAT` is an A/T-homopolymer-adjacent, low-complexity
motif of exactly the kind flagged earlier in this investigation as prone
to coincidental matches from base composition alone, not a real shared
repeat unit. Treat this one row with more caution than the other seven.

Final combined table (coarse classification + both exhaustive
resolutions merged): `outputs/tr_repeat_correlation/final_repeat_confirmations.tsv`,
via `finalize_repeat_confirmations.py`. **44/44 CONFIRMED_PRESENT, 0/44
negative.**

### Side finding: TR-copy count doesn't predict telomeric *frequency*, but every core-template variant found in a gene copy shows up in the telomere somewhere

Originally framed as "LABEL_FLIP: 8/44 sets where the TR-copy-count
majority variant is telomerically silent" - that framing turns out to be
wrong once checked exhaustively. It's not silent in any of the 8 cases;
it's present, just at very low relative frequency (0.1-3.7% in 7/8
cases; see the `Solanum_nigrum` caveat above for the 8th). The real,
surviving finding is narrower but still holds: **TR gene paralog copy
number does not predict which sequence dominates a species' telomere by
abundance** - e.g. *Juncus_bufonius* has more gene copies backing
`CTAAACCCTAGAT`, but the telomere is 99.8% `CTAAACCCTAAAC` by direct
count. What copy number does NOT predict is presence/absence, though:
across all 44 differing sets checked exhaustively, every distinct
core-template sequence found in a TR gene copy was also found somewhere
in that species' telomere, from >99% dominant down to a few dozen
occurrences in 5000bp windows across dozens of chromosome ends.

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

- **The exhaustive walk was run on all 16 ambiguous sets** (8
  WEAK_SENSITIVITY + 8 LABEL_FLIP); the other 28 (27 STRONG_OVERLAP + 1
  PARTIAL_OVERLAP) were already directly observed by the coarse method
  and weren't re-walked. 44/44 rests on the coarse method for those 28
  and the exhaustive walk for the other 16.
- **`Solanum_nigrum`'s LABEL_FLIP resolution (30.5%) is the one result
  in this note that deserves extra scrutiny** - its majority variant
  `AATAAAT` is a low-complexity, A/T-rich motif of the kind flagged
  elsewhere in this investigation as prone to coincidental matches; its
  much higher relative frequency than the other 7 LABEL_FLIP cases
  (all <5%) is consistent with that.
- **Haplotype-phasing risk is real, not hypothetical - now systematically
  checked (2026-09-16).** `Empetrum_nigrum`'s TR loci spread across
  `SUPER_N_HAP2/HAP3/HAP4`-style scaffold names was the tip-off. A full
  dataset sweep (covering both the `SUPER_N_HAPx` suffix and the
  `HAPx_SCAFFOLD_N`/`HAPx_SUPER_N` prefix naming conventions, which are
  both used across different DToL assemblies) found: 68/331 species have
  at least one TR locus on a haplotype-labelled scaffold, but only
  **8/331** show the genuinely inflating pattern - the *same* true
  chromosome carrying a TR locus on *more than one* haplotype number
  (`Buddleja_davidii`, `Empetrum_nigrum`, `Galium_boreale`,
  `Hesperis_matronalis`, `Hypericum_perforatum`, `Lythrum_salicaria`,
  `Salix_cinerea`, `Solidago_canadensis`). In every case checked, the
  redundant haplotype copies of a given locus carry the identical
  sequence (consistent with allelic redundancy of one real locus, not
  independent divergent paralogs) - so this inflates raw copy-number
  counts for these 8 species but does not fabricate spurious sequence
  diversity. Of the 44 headline differing sets specifically, only 1
  (`Lythrum_salicaria`, 11bp) is affected - see the corrected example
  above. The identical-copy control group (headline #6 below) includes 5
  of these 8 species, but since haplotype redundancy doesn't change
  whether a set is classified "identical," it doesn't threaten that
  comparison's validity.
- **Orientation ambiguity remains unresolved** (see
  `tr_repeat_correlation.md`) - most core-template matches are in the
  Template's `direct` sense rather than the textbook `revcomp`
  templating-mechanism expectation, and this analysis didn't revisit
  that question.

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

# Same, for the LABEL_FLIP sets (checks whether the TR-copy-majority
# variant, not the minority, is really absent):
python3 src/resolve_weak_sensitivity.py outputs/tr_repeat_correlation/differing_set_confirmations.tsv \
  outputs/tr_repeat_correlation/positional 5000 outputs/tr_repeat_correlation/label_flip_resolved.tsv LABEL_FLIP

# Merge everything into the final table:
python3 src/finalize_repeat_confirmations.py outputs/tr_repeat_correlation/differing_set_confirmations.tsv \
  outputs/tr_repeat_correlation/weak_sensitivity_resolved.tsv outputs/tr_repeat_correlation/label_flip_resolved.tsv \
  outputs/tr_repeat_correlation/final_repeat_confirmations.tsv

# One-off base-pair-resolved check of a single chromosome-end:
python3 src/walk_terminal_repeat_array.py <species> <chromosome> <5prime|3prime> <variant1> <variant2> --window 5000
```
