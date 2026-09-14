# Structure of TR core-template variation: per-species breakdown

Follow-on to `tr_core_template_variation_and_telomeres.md`, which established
that all 44 differing (species, core-length) sets have both their variants
genuinely present in the telomere. This note asks the finer-grained
question directly: **for each of the 33 species, how many TR gene copies
are there, which are variable, where exactly do they vary, which sequence
actually dominates the real telomere array, and - where copies are
variable - is the minor variant scattered randomly through the array or
does it form a block/HOR-like structure?**

## Headline

- **41/44 (93%) sets show clustered, non-random structure** - the
  non-dominant sequence is essentially never a lone, isolated substitution
  surrounded by the majority sequence. It occurs in contiguous runs, from
  small local patches (2-34 units) up to enormous single blocks (up to
  299 contiguous units). Only 1/44 set had a clean isolated single
  substitution with enough data to call it "random"; 2/44 had only a
  single minor-variant occurrence genome-wide (too little data for any
  verdict).
- **25/44 (57%) sets have a telomere-dominant sequence that differs from
  the TR-copy-count dominant** - much higher than the 8/44 ("LABEL_FLIP")
  found by the earlier top-2-only analysis, because this pass tiles with
  *every* distinct variant at a locus, not just the top 2, and catches
  cases where a third/fourth minor TR-copy variant turns out to be the
  real telomeric majority (see `Ajuga_chamaepitys` core-length-7 example
  below).
- **A genuinely novel structural class turned up: regular, short-period
  alternation between two variants** in ~5 sets (e.g. `Centaurium_intermedium`:
  a clean `...AACCCTAA-CCTAAACC-AACCCTAA-CCTAAACC...` alternation over
  dozens of units) - not random, not one big block either, but a tight,
  repeating two-unit period. This is the closest thing in this dataset to
  a literal higher-order-repeat (HOR) signature.
- **Important caveat discovered along the way**: several "differing"
  cores are periodic telomere-repeat-family motifs, so raw (non-rotation-
  corrected) sequence comparison overstates how different two cores look.
  6/55 non-modal variants are *pure phase reads of the identical sequence*
  (0 true substitutions once rotation is accounted for) - see below.

## Data and new scripts

Three new scripts, run over the 44 differing sets from
`core_template_sets.tsv` / `differing_set_confirmations.tsv`:

| script | role |
|---|---|
| `tr_core_template_variant_positions.py` | full (not just top-2) per-set variant breakdown from `all_species.tsv`, with rotation-corrected within-core diff positions |
| `tr_core_template_array_structure.py` | exhaustive walk (every main chromosome, both ends - same machinery as `resolve_weak_sensitivity.py`) tiling with *every* distinct variant per set; reports true telomere-wide dominance and block/runs-test structure of the non-dominant sequence |
| `tr_core_template_array_structure_summary.py` | categorizes each set by minor-variant frequency band and structure verdict |

Outputs: `outputs/tr_repeat_correlation/core_template_variant_positions.tsv`,
`telomere_array_structure.tsv`, `telomere_array_structure_summary.tsv`,
and a combined `tr_core_template_array_structure_master.tsv` used to build
the table below.

## Caveat first: periodic motifs need rotation-corrected comparison

Most of these cores are built from a short (~3bp) repeating unit (the
`AAC`/`CCT`-family), so two chunks compared in their raw, as-extracted
reading frame can look far more different than they really are - e.g.
`AACCCT` vs `CCCTAA` is a naive 5/6-position mismatch, but `CCCTAA` is
exactly `AACCCT` read starting 2 bases later (0 true substitutions).
Searching all cyclic rotations of each non-modal variant against the
fixed modal frame and keeping the minimum-Hamming-distance alignment
drops the typical apparent difference from 1-13 raw mismatches down to
1-2 true substitutions in the great majority of cases (43/55 non-modal
variants improve), and reveals that **6/55 are pure phase reads of the
identical sequence - not point-mutant paralogs at all**:

| species | core length | modal | "variant" | true differences |
|---|---|---|---|---|
| Anagallis_arvensis | 6 | AACCCT | CCCTAA | 0 (phase shift only) |
| Platanus_x_hispanica | 7 | CCTAAAC | AAACCCT | 0 (phase shift only) |
| Potentilla_anglica | 12 | ACTCTAAACTCT | CTCTAAACTCTA | 0 (phase shift only) |
| Acaena_novae_zelandiae | 14 | AAACCCTAAACCCT | AACCCTAAACCCTA | 0 (phase shift only) |
| Acaena_ovalifolia | 14 | AAACCCTAAACCCT | AACCCTAAACCCTA | 0 (phase shift only) |
| Platanus_x_hispanica | 6 | AAACCT | CTAAAC | 0 (phase shift only) |

For 4 of these 6, this is the *only* non-modal variant in the set -
`Platanus_x_hispanica` (both its sets), `Acaena_novae_zelandiae`, and
`Acaena_ovalifolia` are marked `pure_rotation_only` below, and their
"structure" results (large blocks, near-50/50 splits) most likely
reflect **chromosome-to-chromosome reading-frame/registration
differences of one homogeneous repeat**, not two competing sequences -
confirmed directly for `Platanus_x_hispanica` core-length-7: a hand
inspection of `SUPER_20:5prime` found the array reads 100% in the
`CCTAAAC` phase for its full length, never switching to `AAACCCT` -
whichever phase a chromosome's array happens to start in, it stays
locked there (the greedy tiler advances exactly one repeat-unit length
per match, so it can't drift once registered); different chromosomes
apparently commit to different starting phases, which is what produces
the genome-wide near-50/50 split and huge single blocks. This does **not**
affect the earlier 44/44 telomere-confirmation headline (that's about
real, literal DNA content, true regardless of framing) - it only means
"variant" isn't the right word for these 4 sets specifically.

This does not undermine the genuine point-mutant cases (the other 49/55
non-modal variants do carry 1-5 true substitutions after rotation
correction) - it just means positions should be read from the
rotation-corrected column, not the raw one.

## Per-species/per-set table

Sorted by species. `n_loci`/`n_distinct` = TR gene copies and distinct
core sequences at that core length (from `core_template_sets.tsv`).
`true_diff` = rotation-corrected substitution count of the top minority
variant vs the modal (0 = pure phase rotation, see above).
`telomere_dominant` = the sequence actually most common in the real
telomere, from the exhaustive walk (not necessarily the TR-copy modal).
`minor_freq` = fraction of all telomere-matched repeat units belonging to
non-dominant sequence(s). `structure` = block/runs-test verdict.

| species | core len | n_loci | n_distinct | true_diff | TR-copy modal | telomere-dominant | copy≠telomere dom? | minor freq | band | structure |
|---|---|---|---|---|---|---|---|---|---|---|
| Acaena_novae_zelandiae | 14 | 2 | 2 | 0* | AAACCCTAAACCCT | AACCCTAAACCCTA | yes | 44.7% | CO-DOMINANT | block of 106 |
| Acaena_ovalifolia | 14 | 2 | 2 | 0* | AAACCCTAAACCCT | AACCCTAAACCCTA | yes | 43.9% | CO-DOMINANT | block of 129 |
| Ajuga_chamaepitys | 9 | 7 | 2 | 2 | AAACTAAAC | AAACTAAAC | no | 3.9% | RARE | block of 3 |
| Ajuga_chamaepitys | 7 | 6 | 4 | 3 | CTAAACT | **CTAAACC** (=canonical repeat, phase-shifted) | **yes** | 6.7% | RARE | block of 34 |
| Ajuga_chamaepitys | 6 | 3 | 2 | 1 | AACCCT | AACCCT | no | 0.2% | RARE | block of 2 |
| Ajuga_chamaepitys | 8 | 2 | 2 | 2 | CTAAACCT | CTAAACCT | no | 38.0% | CO-DOMINANT | block of 33 |
| Ajuga_chamaepitys | 12 | 2 | 2 | 3 | AAAACCCTAAAC | AAAACCCTAAAC | no | 45.6% | CO-DOMINANT | block of 29 |
| Ajuga_reptans | 11 | 2 | 2 | 2 | ACCCTAAACCT | AACCCTGAACC | yes | 48.1% | CO-DOMINANT | block of 37 |
| Alchemilla_mollis | 12 | 2 | 2 | 2 | TAAACCCTAAAC | AACCCTAAACCC | yes | 47.0% | CO-DOMINANT | block of 106 |
| Anagallis_arvensis | 6 | 17 | 6 | 0* | AACCCT | CCCTAA | yes | 53.5% | CO-DOMINANT | block of 77 |
| Anagallis_arvensis | 9 | 4 | 2 | 2 | CCCTAAACC | AAACCCTAA | yes | 29.9% | MINORITY | block of 4 |
| Anagallis_arvensis | 11 | 2 | 2 | 2 | AACCCTAAACC | AACCCTAAACC | no | 45.3% | CO-DOMINANT | block of 27 |
| Arctium_lappa | 13 | 4 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 2.2% | RARE | block of 3 |
| Arctium_minus | 13 | 3 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 1.3% | RARE | isolated (no blocks) |
| Carpinus_betulus | 13 | 5 | 2 | 2 | TAAACCCTAAAAG | CTAAACCCTAAAC | yes | 1.4% | RARE | block of 2 |
| Centaurium_intermedium | 8 | 2 | 2 | 1 | CCTAAACC | AACCCTAA | yes | 49.9% | CO-DOMINANT | **regular alternation** (block cap 3, z=+21) |
| Cornus_sanguinea | 13 | 2 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 1.5% | RARE | block of 4 |
| Crataegus_laevigata | 6 | 2 | 2 | 1 | CAACCT | TAAACC | yes | 0.5% | RARE | block of 4 |
| Glaucium_flavum | 14 | 2 | 2 | 1 | GAAAACCCTACCCG | GAAAACCCTACCCG | no | 3.3% | RARE | block of 5 |
| Hedera_helix | 9 | 2 | 2 | 2 | TAAACCCTA | CCTAAACCC | yes | 49.4% | CO-DOMINANT | **regular alternation** (block cap 6, z=+22) |
| Helianthemum_oelandicum | 11 | 5 | 3 | 2 | ACCCCTAACCC | ACCCTAACCCT | yes | 52.5% | CO-DOMINANT | block of 26 |
| Inula_conyza | 13 | 4 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 2.4% | RARE | block of 4 |
| Juncus_bufonius | 13 | 3 | 2 | 2 | CTAAACCCTAGAT | CTAAACCCTAAAC | yes | 0.2% | RARE | block of 2 |
| Linaria_repens | 9 | 4 | 2 | 1 | CTAAACCCT | CCCTAAACC | yes | 49.7% | CO-DOMINANT | **regular alternation** (block cap 4, z=+21) |
| Lycopus_europaeus | 12 | 2 | 2 | 2 | CTAAACCCTACC | AAACCCTAAACC | yes | 0.1% | RARE | 1 occurrence only |
| Lythrum_hyssopifolia | 10 | 2 | 2 | 1 | AACCCTAAAC | ACCCTAAACC | yes | 48.1% | CO-DOMINANT | block of 148 |
| Lythrum_portula | 12 | 2 | 2 | 2 | AACCCTAAACCC | AACCCTAACCCT | yes | 32.2% | MINORITY | block of 9 |
| Lythrum_salicaria | 11 | 6 | 3 | 1 | AACCCTAAACC | ACCCTAAACCC | yes | 45.5% | CO-DOMINANT | block of 238 |
| Parietaria_judaica | 11 | 2 | 2 | 1 | ACCCTAAACCC | AACCCTAAACC | yes | 46.2% | CO-DOMINANT | block of 225 |
| Platanus_x_hispanica | 7 | 6 | 2 | 0* | CCTAAAC | CCTAAAC | no | 49.0% | CO-DOMINANT | block of 299 |
| Platanus_x_hispanica | 6 | 2 | 2 | 0* | AAACCT | CTAAAC | yes | 0.9% | RARE | block of 8 |
| Pontederia_cordata | 15 | 2 | 2 | 1 | TAAACCCTAAACCCT | AAACCCTAAACCCTA | yes | 48.0% | CO-DOMINANT | **regular alternation** (block cap 3, z=+16) |
| Potentilla_anglica | 13 | 6 | 2 | 3 | ACTCTAAACTCTA | ACTCTAAACTCTA | no | 28.9% | MINORITY | block of 16 |
| Potentilla_anglica | 12 | 3 | 3 | 0* | ACTCTAAACTCT | TAAACCCTAAAC | yes | 62.0% | CO-DOMINANT | block of 41 |
| Potentilla_indica | 14 | 3 | 3 | 1 | TAAACCCTAAACCT | AAACCCTAAACCCT | yes | 48.0% | CO-DOMINANT | block of 47 |
| Potentilla_indica | 11 | 2 | 2 | 1 | CTAAACTCTAA | CTAAACTCTAA | no | 30.0% | MINORITY | block of 12 |
| Primula_scotica | 8 | 2 | 2 | 1 | AACCCTAA | AACCCTAA | no | 49.8% | CO-DOMINANT | **regular alternation** (block cap 5, z=+19.5) |
| Ranunculus_repens | 12 | 5 | 2 | 2 | TGAACCCTGAAC | TGAACCCTGAAC | no | 12.4% | MINORITY | block of 7 |
| Solanum_nigrum | 11 | 12 | 3 | 3 | CTGAACCCTGA | CTGAACCCTGA | no | 56.1% | CO-DOMINANT | block of 268 |
| Solanum_nigrum | 7 | 2 | 2 | 4 | AATAAAT | ACCTGAA | yes | 30.5% | MINORITY | block of 12 (low-complexity motif - caveat) |
| Solanum_nigrum | 12 | 2 | 2 | 2 | TAAACCCTAAAC | TAAACCCTAAAC | no | 24.9% | MINORITY | block of 23 |
| Tellima_grandiflora | 11 | 3 | 2 | 2 | AACCCTAAACC | AACCCTAAACC | no | 4.8% | RARE | block of 22 |
| Tilia_cordata | 13 | 2 | 2 | 1 | CAAACCCTAAACC | TAAACCCTAAACC | yes | 1.3% | RARE | block of 10 |
| Trocdaris_verticillatum | 9 | 5 | 2 | 5 | ACCCTAACC | ACCCTAACC | no | 0.7% | RARE | 1 occurrence only |

`*` = pure phase rotation (0 true substitutions, see caveat above).

## Which is dominant in the telomere?

25/44 sets (57%) have a telomere-dominant sequence that is *not* the
TR-copy-count dominant. TR gene paralog count simply does not predict
which sequence a species' telomere array actually runs on. The most
striking example: `Ajuga_chamaepitys` core-length-7 - the TR-copy modal
is `CTAAACT` (3/6 gene copies), but the real telomere is 93% `CTAAACC`,
which is exactly the canonical plant telomere repeat (`TTTAGGG`'s
reverse complement, `CCCTAAA`, read in a different phase) - a variant
that was carried by only 1/6 gene copies and was invisible to the
original top-2-only pipeline entirely.

## Where do they vary within the core?

Rotation-corrected substitution counts for the 49 genuine (non-pure-
rotation) non-modal variants: mostly 1-2 true differences (41/49, 84%),
occasionally up to 5. No indels within any set (by construction - the
comparable-set grouping buckets loci by exact core length, so this isn't
an independent finding, just confirms every within-set comparison really
is length-matched). A fractional-position histogram (substitution position /
core length) shows more mismatches concentrated near the two ends of the
core (positions near 0.0 and 1.0) than the middle - but this is reported
with a caveat, not as a clean biological claim: the ends of an extracted
"chunk" are exactly where its own boundary-detection algorithm
(`tr_repeat_correlation.py`'s best-match search) is least well anchored,
so some of this apparent edge-enrichment could be a boundary-detection
artifact rather than the true templating region genuinely being more
tolerant of substitution at its edges. Distinguishing the two would need
an independent boundary check (e.g. cross-referencing against the
telomere-independent Template-boundary prediction in
`tr_template_boundary_prediction.md`), not done here.

## What is the array structure - random, clustered, or a pattern?

Three regimes, not one:

**RARE + blocky (16/44 sets, minor freq <10%).** The classic picture -
one sequence overwhelmingly dominant, a genuine point-mutant variant
present at low frequency. But it is essentially never isolated: 13/16
form contiguous multi-unit blocks (2-34 units) somewhere in the genome,
not scattered single substitutions. Only `Arctium_minus` showed a clean
isolated pattern with enough data to call it consistent with random
placement; `Lycopus_europaeus` and `Trocdaris_verticillatum` had only 1
occurrence genome-wide (no verdict possible). This argues against simple
recurrent point mutation as the source of the minor variant, and for
something like a localized tandem duplication/expansion (or gene
conversion copying a short stretch) having produced a multi-unit patch
at one place in the array, which then persists.

**MINORITY/CO-DOMINANT + large blocks (28/44 sets, minor freq 12-62%).**
Both sequences are substantial fractions of the array, arranged in large
contiguous domains (blocks of dozens to hundreds of units) rather than
finely interspersed. For 3 of the 4 pure-phase-rotation sets
(`Acaena_novae_zelandiae`, `Acaena_ovalifolia`, `Platanus_x_hispanica`
core-length-7 - the 4th, `Platanus_x_hispanica` core-length-6, falls in
the RARE bucket above) this most likely reflects between-chromosome
registration differences of one repeat (see caveat above), not two
competing sequences - but the large majority of this bucket (25/28) are
genuine substitution pairs, where large single-variant domains of
hundreds of contiguous repeat units is itself a striking,
previously-undocumented structural finding.

**Regular short-period alternation (5/44 sets: `Centaurium_intermedium`,
`Hedera_helix`, `Linaria_repens`, `Pontederia_cordata`, `Primula_scotica`).**
All near-50/50 splits, but with small block caps (3-6 units) *and* a
strongly positive runs-test z-score (+16 to +22) - far more alternation
between the two sequences than random interspersion would produce. This
is the closest thing in the dataset to a literal composite/HOR unit.
Direct inspection of `Centaurium_intermedium` (`SUPER_1_HAP1:3prime`)
shows a textbook example:

```
...AACCCTAA / CCTAAACC / AACCCTAA / CCTAAACC / AACCCTAA / CCTAAACC...
```

18 consecutive units alternating almost perfectly, repeated at multiple
independent chromosome ends. Whether this reflects a true two-variant
compound repeat unit encoded that way in the genome, versus some
artifact of how a single physical unit maps to two candidate strings,
isn't resolved here - worth a closer, sequence-level look before
treating it as a confirmed HOR.

## Caveats

- **Structural verdicts rest on a single genome-wide "best" chromosome-
  end for the formal runs-test z-score** (whichever has the most
  minor-variant hits) - not pooled across all chromosome-ends
  statistically. `max_minor_block`, by contrast, IS a genuine genome-wide
  maximum (checked across every chromosome, both ends).
- **Low-count cases have no statistical power.** 2/44 sets had only 1
  minor-variant occurrence genome-wide; their "structure" is unknowable,
  not confirmed random.
- **`Solanum_nigrum`'s core-length-7 set (`AATAAAT`/`ACCTGAA`) involves a
  low-complexity, A/T-rich motif** flagged earlier in this investigation
  as prone to coincidental matches - treat its 30.5% minor frequency and
  block structure with extra caution.
- **The 4 pure-phase-rotation sets' "structure" likely reflects
  chromosome-to-chromosome reading-frame registration of one repeat, not
  two competing sequences** - see the caveat section above.
- **The greedy left-to-right tiler stays phase-locked once it starts
  matching** (it advances by exactly one repeat-unit length per match),
  so within one contiguous, correctly-phased stretch it cannot itself
  register-shift - a real single-base indel in the underlying sequence
  would break the tiling into a gap until the next chance re-sync. This
  is a general limitation of the method (already noted for tidk vs. this
  approach in `tr_repeat_correlation.md`), not specific to this analysis.
- **Position-within-core histogram's edge-enrichment may be partly a
  chunk-boundary-detection artifact**, not a clean biological signal -
  see "Where do they vary" above.

## Reproducing

```bash
# Full (not top-2-only) variant breakdown + rotation-corrected diff positions:
python3 src/tr_core_template_variant_positions.py outputs/tr_repeat_correlation/all_species.tsv \
  outputs/tr_repeat_correlation/core_template_variant_positions.tsv

# Exhaustive telomere-wide dominance + block/runs-test structure (every main
# chromosome, both ends, one genome decompression per species):
python3 src/tr_core_template_array_structure.py outputs/tr_repeat_correlation/core_template_variant_positions.tsv \
  outputs/tr_repeat_correlation/positional 5000 outputs/tr_repeat_correlation/telomere_array_structure.tsv

# Categorize into frequency bands + structure verdicts:
python3 src/tr_core_template_array_structure_summary.py outputs/tr_repeat_correlation/telomere_array_structure.tsv \
  outputs/tr_repeat_correlation/telomere_array_structure_summary.tsv

# One-off illustrative render of a specific chromosome-end:
python3 src/walk_terminal_repeat_array.py Centaurium_intermedium SUPER_1_HAP1 3prime CCTAAACC AACCCTAA --window 300
```
