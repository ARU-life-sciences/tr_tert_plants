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

- **A fixed 5000bp terminal window can include subtelomeric/interstitial
  sequence, not just the true telomere** - plant telomere length varies
  ~0.3-200kb+ across species (Watson & Riha 2010 FEBS Lett), so a fixed
  window could sit entirely within the telomere for a long-telomere
  species or well past it for a short one. All structure verdicts below
  are computed **restricted to "telomere-proper"** - the contiguous run
  of matches starting exactly at the true chromosome terminus, up to the
  first >50bp break - not the raw full window. This is a real correction,
  not a refinement: 3/44 differing sets and 17/62 control sets (see below)
  lose their apparent minor-variant signal entirely once restricted,
  meaning it was subtelomeric noise, not real telomere content.
- **35/44 (80%) sets that retain a confirmed minor variant show
  non-random structure** (29 blocky/clustered + 6 regular-alternation);
  only 6/44 show isolated substitutions consistent with random placement,
  and 3/44 show no minor-variant presence at all once restricted to
  telomere-proper.
- **Controlled comparison, the main new result of this note**: does this
  structure need TR paralog divergence, or would it show up anyway
  (e.g. from DNA-level recombination independent of any TR gene, per
  Belyayev et al. 2023)? Tested directly by running the identical
  exhaustive-walk method on 62 sets from species with **fully identical**
  TR gene copies (no paralog divergence at all), using an independently
  tidk-discovered second candidate (blind to TR content) instead of a
  real paralog. Result:

  | | TR-paralog-differing (44) | TR-identical control (62) | Fisher's exact p |
  |---|---|---|---|
  | Regular alternation (candidate HOR) | 6/44 (13.6%) | **0/62** | **0.004** |
  | Any confirmed minor variant | 41/44 (93.2%) | 45/62 (72.6%) | 0.011 |
  | Non-random structure, given present | 35/41 (85.4%) | 30/45 (66.7%) | 0.049 |

  Ordinary "blocky" clustering happens at a real baseline rate even
  without any TR paralog divergence (consistent with a general DNA-level
  mechanism, as Belyayev et al. propose) - but **regular short-period
  alternation was found only in TR-paralog-differing loci, in none of the
  62 controls**. This is a between-species comparison (the 6
  regular-alternation species have no matched identical-copy locus that
  passed into the control group), not a same-species paired design.
- **25/44 (57%) sets have a telomere-dominant sequence that differs from
  the TR-copy-count dominant** - much higher than the 8/44 ("LABEL_FLIP")
  found by the earlier top-2-only analysis, because this pass tiles with
  *every* distinct variant at a locus, not just the top 2.
- **Important caveat discovered along the way**: several "differing"
  cores are periodic telomere-repeat-family motifs, so raw (non-rotation-
  corrected) sequence comparison overstates how different two cores look.
  6/55 non-modal variants are *pure phase reads of the identical sequence*
  (0 true substitutions once rotation is accounted for) - see below.

## Data and new scripts

| script | role |
|---|---|
| `tr_core_template_variant_positions.py` | full (not just top-2) per-set variant breakdown from `all_species.tsv`, with rotation-corrected within-core diff positions |
| `tr_core_template_array_structure.py` | exhaustive walk (every main chromosome, both ends - same machinery as `resolve_weak_sensitivity.py`) tiling with *every* distinct variant per set over the 44 differing sets; reports telomere-wide dominance and block/runs-test structure, BOTH over the full extraction window and restricted to telomere-proper (see below) |
| `tr_core_template_array_structure_summary.py` | categorizes each set (from the telomere-proper columns) by minor-variant frequency band and structure verdict, with the z-test-first classification fix (see below) |
| `find_identical_set_second_candidates.py` | for species with fully identical TR copies, finds a plausible independent second repeat candidate from that species' own tidk-explore output (blind to TR content) - the control group's input |
| `tr_identical_set_array_structure.py` | same exhaustive-walk machinery applied to the identical-copy control group; defines "main chromosome" directly from genome assembly length (`samtools faidx`) rather than reusing a tidk-search window count, since there's no pre-existing scan to reuse for these species |

Outputs: `outputs/tr_repeat_correlation/core_template_variant_positions.tsv`,
`telomere_array_structure.tsv` (44 differing sets, full_*/proper_* columns),
`telomere_array_structure_summary.tsv`, `identical_set_second_candidates.tsv`,
`identical_set_variant_positions.tsv`, `identical_set_array_structure.tsv`
(62 control sets), `identical_set_array_structure_summary.tsv`, and a
combined `tr_core_template_array_structure_master.tsv` used to build the
table below.

## Methodology corrections made while building the control-group comparison

**Telomere-proper boundary restriction.** The original pass tiled the
full fixed 5000bp terminal extraction window. Plant telomere length
varies enormously across species (documented range ~0.3kb to 200kb+ -
Watson & Riha 2010 FEBS Lett; Procházková Schrumpfová et al. 2016 Front
Plant Sci) - too wide a range to assume the window is either safely
inside or safely past the true telomere for any given species. Belyayev
et al. 2023 also note interstitial telomeric repeats (ITRs) are *more*
degenerate than terminal arrays, so pooling whatever lies beyond the true
telomere risks mistaking ITR noise for real terminal structure. Fixed by
`telomere_proper_prefix()`: walk from the true chromosome terminus and
stop at the first gap >50bp in matches - adaptive to whatever the real
array length is (a long, dense telomere is never truncated; a short one
is), not tied to an assumed length. Real effect, not a formality: 3/44
differing sets and 17/62 control sets lose their entire minor-variant
signal once restricted - it was subtelomeric noise.

**Verdict-classification bug.** The original per-set verdict checked "is
there any block >=2 units" *before* checking the formal runs-test
statistic. A strongly regular, short-period alternating pattern (many
short blocks, switching back and forth far more than random - a strongly
*positive* z) technically also satisfies "a block of >=2 exists
somewhere," so it was silently mislabelled "CLUSTERED" instead of the
more informative "REGULAR_ALTERNATION" - exactly the pattern most
relevant to a genuine HOR claim. Fixed by checking the z-test first:
significant positive z -> REGULAR_ALTERNATION, significant negative z or
(when the test is unpowered) a directly observed block -> CLUSTERED.
This reclassification is why the regular-alternation count moved from an
earlier informal "~5" (found by eyeballing high z-scores) to a properly
threshold-derived 6/44 (z > 1.96), adding `Anagallis_arvensis`
core-length-9 (z=+3.6, weaker but still significant) to the original 5.

**Identical-copy control group.** Built from the 223 comparable sets
where all TR gene copies share an identical core (`n_distinct==1` in
`core_template_sets.tsv`) - restricted to core lengths 6-15bp for
comparability, giving 215 candidate sets / 194 species. For each,
`find_identical_set_second_candidates.py` searches that species' own
`outputs/tidk/<species>.tidk.tsv` (independent tandem-repeat discovery,
blind to TR gene content) for a same-length, rotation-corrected
low-edit-distance (<=3bp), non-trivial-count second candidate - 62 sets
(60 species) had one. The same exhaustive walk and telomere-proper
restriction was then applied to each (species, dominant, candidate)
triple. Two things distinguish this from simply reusing the
differing-set pipeline: (1) "main chromosome" is defined directly from
genome assembly length via `samtools faidx`, not a tidk-search window
count, since there was no pre-existing scan to reuse for these species
(arguably more principled than the differing-set pipeline's approach, not
just a shortcut - both target the same concept of "large, well-assembled
scaffold"); (2) candidates are inherently noisier a priori than the
differing sets' TR-gene-verified variants, since they come from an
unconstrained same-length + small-edit-distance search over each
species' full tidk candidate list, not from an independently-validated
biological paralog - the decisive filter is real telomere-proper presence
in the walk itself (17/62 candidates failed this and were excluded from
the final comparison), matching the "always confirm by walking" principle
used throughout this whole investigation.

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
| Acaena_novae_zelandiae | 14 | 2 | 2 | 0* | AAACCCTAAACCCT | AACCCTAAACCCTA | yes | 37.5% | CO-DOMINANT/BALANCED | CLUSTERED (runs test, block 106, z=-6.0) |
| Acaena_ovalifolia | 14 | 2 | 2 | 0* | AAACCCTAAACCCT | AACCCTAAACCCTA | yes | 43.9% | CO-DOMINANT/BALANCED | BLOCKY (longest run 129 units) |
| Ajuga_chamaepitys | 9 | 7 | 2 | 2 | AAACTAAAC | AAACTAAAC | no | 0.0% | ABSENT | ABSENT |
| Ajuga_chamaepitys | 7 | 6 | 4 | 3 | CTAAACT | CTAAACC | yes | 7.0% | RARE | CLUSTERED (runs test, block 34, z=-13.1) |
| Ajuga_chamaepitys | 6 | 3 | 2 | 1 | AACCCT | AACCCT | no | 0.2% | RARE | ISOLATED (never adjacent to itself) |
| Ajuga_chamaepitys | 8 | 2 | 2 | 2 | CTAAACCT | CTAAACCT | no | 29.0% | MINORITY | CLUSTERED (runs test, block 6, z=-2.5) |
| Ajuga_chamaepitys | 12 | 2 | 2 | 3 | AAAACCCTAAAC | AAAACCCTAAAC | no | 42.9% | CO-DOMINANT/BALANCED | BLOCKY (longest run 3 units) |
| Ajuga_reptans | 11 | 2 | 2 | 2 | ACCCTAAACCT | ACCCTAAACCT | no | 20.0% | MINORITY | ISOLATED (never adjacent to itself) |
| Alchemilla_mollis | 12 | 2 | 2 | 2 | TAAACCCTAAAC | AACCCTAAACCC | yes | 47.7% | CO-DOMINANT/BALANCED | BLOCKY (longest run 106 units) |
| Anagallis_arvensis | 6 | 17 | 6 | 0 | AACCCT | CCCTAA | yes | 50.1% | CO-DOMINANT/BALANCED | CLUSTERED (runs test, block 77, z=-13.6) |
| Anagallis_arvensis | 9 | 4 | 2 | 2 | CCCTAAACC | AAACCCTAA | yes | 30.9% | MINORITY | REGULAR_ALTERNATION (block cap 4, z=+3.6) |
| Anagallis_arvensis | 11 | 2 | 2 | 2 | AACCCTAAACC | CCCTAAACCCT | yes | 44.4% | CO-DOMINANT/BALANCED | BLOCKY (longest run 17 units) |
| Arctium_lappa | 13 | 4 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 1.2% | RARE | BLOCKY (longest run 2 units) |
| Arctium_minus | 13 | 3 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 0.3% | RARE | ISOLATED (never adjacent to itself) |
| Carpinus_betulus | 13 | 5 | 2 | 2 | TAAACCCTAAAAG | CTAAACCCTAAAC | yes | 0.1% | RARE | ISOLATED (never adjacent to itself) |
| Centaurium_intermedium | 8 | 2 | 2 | 1 | CCTAAACC | CCTAAACC | no | 49.9% | CO-DOMINANT/BALANCED | REGULAR_ALTERNATION (block cap 4, z=+21.3) |
| Cornus_sanguinea | 13 | 2 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 1.4% | RARE | CLUSTERED (runs test, block 2, z=-2.2) |
| Crataegus_laevigata | 6 | 2 | 2 | 1 | CAACCT | TAAACC | yes | 0.0% | RARE | CLUSTERED (runs test, block 2, z=-11.5) |
| Glaucium_flavum | 14 | 2 | 2 | 1 | GAAAACCCTACCCG | GAAAACCCTACCCG | no | 1.0% | RARE | CLUSTERED (runs test, block 2, z=-3.1) |
| Hedera_helix | 9 | 2 | 2 | 2 | TAAACCCTA | CCTAAACCC | yes | 49.6% | CO-DOMINANT/BALANCED | REGULAR_ALTERNATION (block cap 4, z=+21.7) |
| Helianthemum_oelandicum | 11 | 5 | 3 | 2 | ACCCCTAACCC | ACCCTAACCCT | yes | 52.6% | CO-DOMINANT/BALANCED | CLUSTERED (runs test, block 26, z=-7.9) |
| Inula_conyza | 13 | 4 | 2 | 1 | TAAACCCTAAACC | TAAACCCTAAACC | no | 0.8% | RARE | CLUSTERED (runs test, block 3, z=-5.3) |
| Juncus_bufonius | 13 | 3 | 2 | 2 | CTAAACCCTAGAT | CTAAACCCTAAAC | yes | 0.1% | RARE | ISOLATED (never adjacent to itself) |
| Linaria_repens | 9 | 4 | 2 | 1 | CTAAACCCT | CTAAACCCT | no | 49.8% | CO-DOMINANT/BALANCED | REGULAR_ALTERNATION (block cap 4, z=+20.9) |
| Lycopus_europaeus | 12 | 2 | 2 | 2 | CTAAACCCTACC | AAACCCTAAACC | yes | 0.1% | RARE | ISOLATED (never adjacent to itself) |
| Lythrum_hyssopifolia | 10 | 2 | 2 | 1 | AACCCTAAAC | ACCCTAAACC | yes | 48.9% | CO-DOMINANT/BALANCED | BLOCKY (longest run 148 units) |
| Lythrum_portula | 12 | 2 | 2 | 2 | AACCCTAAACCC | AACCCTAACCCT | yes | 32.0% | MINORITY | BLOCKY (longest run 9 units) |
| Lythrum_salicaria | 11 | 6 | 3 | 1 | AACCCTAAACC | ACCCTAAACCC | yes | 46.5% | CO-DOMINANT/BALANCED | BLOCKY (longest run 238 units) |
| Parietaria_judaica | 11 | 2 | 2 | 1 | ACCCTAAACCC | AACCCTAAACC | yes | 45.7% | CO-DOMINANT/BALANCED | BLOCKY (longest run 225 units) |
| Platanus_x_hispanica | 7 | 6 | 2 | 0* | CCTAAAC | CCTAAAC | no | 49.3% | CO-DOMINANT/BALANCED | CLUSTERED (runs test, block 299, z=-24.6) |
| Platanus_x_hispanica | 6 | 2 | 2 | 0* | AAACCT | CTAAAC | yes | 0.3% | RARE | CLUSTERED (runs test, block 4, z=-10.8) |
| Pontederia_cordata | 15 | 2 | 2 | 1 | TAAACCCTAAACCCT | AAACCCTAAACCCTA | yes | 49.2% | CO-DOMINANT/BALANCED | REGULAR_ALTERNATION (block cap 3, z=+16.1) |
| Potentilla_anglica | 13 | 6 | 2 | 3 | ACTCTAAACTCTA | ACTCTAAACTCTA | no | 28.5% | MINORITY | CLUSTERED (runs test, block 15, z=-5.9) |
| Potentilla_anglica | 12 | 3 | 3 | 0 | ACTCTAAACTCT | TAAACCCTAAAC | yes | 64.5% | CO-DOMINANT/BALANCED | CLUSTERED (runs test, block 41, z=-4.0) |
| Potentilla_indica | 14 | 3 | 3 | 1 | TAAACCCTAAACCT | TAAACCCTAAACCC | yes | 52.6% | CO-DOMINANT/BALANCED | BLOCKY (longest run 20 units) |
| Potentilla_indica | 11 | 2 | 2 | 1 | CTAAACTCTAA | CTAAACTCTAA | no | 41.8% | CO-DOMINANT/BALANCED | BLOCKY (longest run 12 units) |
| Primula_scotica | 8 | 2 | 2 | 1 | AACCCTAA | CCTAAACC | yes | 50.0% | CO-DOMINANT/BALANCED | REGULAR_ALTERNATION (block cap 13, z=+19.9) |
| Ranunculus_repens | 12 | 5 | 2 | 2 | TGAACCCTGAAC | TGAACCCTGAAC | no | 10.3% | MINORITY | BLOCKY (longest run 7 units) |
| Solanum_nigrum | 11 | 12 | 3 | 3 | CTGAACCCTGA | CTGAACCCTGA | no | 43.4% | CO-DOMINANT/BALANCED | BLOCKY (longest run 35 units) |
| Solanum_nigrum | 7 | 2 | 2 | 4 | AATAAAT | AATAAAT | no | 0.0% | ABSENT | ABSENT |
| Solanum_nigrum | 12 | 2 | 2 | 2 | TAAACCCTAAAC | TAAACCCTAAAC | no | 18.0% | MINORITY | BLOCKY (longest run 3 units) |
| Tellima_grandiflora | 11 | 3 | 2 | 2 | AACCCTAAACC | AACCCTAAACC | no | 0.5% | RARE | BLOCKY (longest run 3 units) |
| Tilia_cordata | 13 | 2 | 2 | 1 | CAAACCCTAAACC | TAAACCCTAAACC | yes | 0.7% | RARE | BLOCKY (longest run 3 units) |
| Trocdaris_verticillatum | 9 | 5 | 2 | 5 | ACCCTAACC | ACCCTAACC | no | 0.0% | ABSENT | ABSENT |

`*` = pure phase rotation (0 true substitutions, see caveat above).

## Which is dominant in the telomere?

25/44 sets (57%) have a telomere-dominant sequence that is *not* the
TR-copy-count dominant (from the telomere-proper-restricted data; the
`copy≠telomere dom?` column above). TR gene paralog count simply does
not predict which sequence a species' telomere array actually runs on.
The most striking example: `Ajuga_chamaepitys` core-length-7 - the
TR-copy modal is `CTAAACT` (3/6 gene copies), but the real telomere is
dominated by `CTAAACC`, which is exactly the canonical plant telomere
repeat (`TTTAGGG`'s reverse complement, `CCCTAAA`, read in a different
phase) - a variant that was carried by only 1/6 gene copies and was
invisible to the original top-2-only pipeline entirely.

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

(All figures below are telomere-proper-restricted; 3/44 sets show no
confirmed minor variant at all and are excluded from the regime counts.)

**RARE + mostly non-random (14/41 sets with confirmed presence, minor
freq <10%).** The classic picture - one sequence overwhelmingly
dominant, a genuine point-mutant variant present at low frequency. But
it is more often than not non-isolated: 9/14 (6 formally CLUSTERED by
the runs test, 3 more show a directly-observed block without enough
data for the formal test) form contiguous multi-unit blocks somewhere in
the genome; 5/14 show a clean isolated pattern consistent with random
placement (`Ajuga_chamaepitys` core-6, `Arctium_minus`, `Carpinus_betulus`,
`Juncus_bufonius`, `Lycopus_europaeus`). This is a real but weaker
signal than the earlier full-window analysis suggested - some of what
looked like clustering in low-frequency cases was subtelomeric noise -
but the balance still favors non-random placement over simple recurrent
point mutation.

**MINORITY/CO-DOMINANT + mostly non-random (27/41 sets, minor freq
10-65%).** Both sequences are substantial fractions of the array; 26/27
show non-random structure (13 blocky, 7 formally clustered, 6 regular
alternation - see below), only 1 isolated. For 3 of the 4 pure-phase-
rotation sets (`Acaena_novae_zelandiae`, `Acaena_ovalifolia`,
`Platanus_x_hispanica` core-length-7 - the 4th, core-length-6, falls in
the RARE bucket) this most likely reflects between-chromosome
registration differences of one repeat (see caveat above), not two
competing sequences - but the large majority (24/27) are genuine
substitution pairs, where large single-variant domains of hundreds of
contiguous repeat units is itself a striking, previously-undocumented
structural finding.

**Regular short-period alternation (6/44 sets, formal z > 1.96:
`Anagallis_arvensis` core-9, `Centaurium_intermedium`, `Hedera_helix`,
`Linaria_repens`, `Pontederia_cordata`, `Primula_scotica`).** All
near-50/50 splits, but with small block caps (3-13 units) *and* a
significantly positive runs-test z-score (+3.6 to +21.7) - far more
alternation between the two sequences than random interspersion would
produce. This is the closest thing in the dataset to a literal
composite/HOR unit, **and the single result in this note found
exclusively among TR-paralog-differing loci - 0/62 in the identical-copy
control group (Fisher's exact p=0.004, see Headline)**. Direct inspection
of `Centaurium_intermedium` (`SUPER_1_HAP1:3prime`) shows a textbook
example:

```
...AACCCTAA / CCTAAACC / AACCCTAA / CCTAAACC / AACCCTAA / CCTAAACC...
```

18 consecutive units alternating almost perfectly, repeated at multiple
independent chromosome ends. Whether this reflects a true two-variant
compound repeat unit encoded that way in the genome, versus some
artifact of how a single physical unit maps to two candidate strings,
isn't fully resolved by sequence inspection alone - but its exclusive
association with TR paralog divergence (never observed in 62
TR-identical controls) is itself strong evidence this is a real,
TR-linked phenomenon rather than a generic mapping artifact, since the
same tiling/mapping method was applied identically to both groups.

## Does this structure need TR paralog divergence, or would it happen anyway?

Belyayev et al. 2023 independently found near-identical block-organized
telomere array structure in *Chenopodium* (unrelated family), via
long-read sequencing and fiber-FISH, and attributed it to DNA-level
recombination between G-/C-rich strands - a mechanism that has nothing
to do with TR gene paralogs (untested in their study). This is a direct,
testable competing hypothesis against this note's framing, and - unlike
Kumawat et al. or Belyayev et al., neither of whom had both TR-divergent
and TR-identical species to compare - our dataset does.

**Control group**: the identical exhaustive-walk method applied to 62
sets from species with fully identical TR gene copies (no paralog
divergence at all), using an independently tidk-discovered second
candidate blind to TR content (`find_identical_set_second_candidates.py`
- see Methodology section above for the selection/validation procedure).

| | TR-paralog-differing (44) | TR-identical control (62) | Fisher's exact p |
|---|---|---|---|
| Regular alternation (candidate HOR) | 6/44 (13.6%) | 0/62 (0%) | 0.004 |
| Any confirmed minor variant | 41/44 (93.2%) | 45/62 (72.6%) | 0.011 |
| Non-random structure, given present | 35/41 (85.4%) | 30/45 (66.7%) | 0.049 |

**Reading**: ordinary blocky/clustered structure happens at a real
baseline rate even without any TR paralog divergence (30/62, 48.4%
overall) - consistent with Belyayev's DNA-recombination mechanism
operating generally across plant telomeres, not requiring TR paralogs.
But regular short-period alternation was found in *zero* of the 62
controls, only in TR-paralog-differing loci. The two mechanisms aren't
mutually exclusive: a general, TR-independent process may explain the
baseline rate of clustering seen in both groups, while TR paralog
divergence specifically adds the fine-grained alternating pattern on top
of that baseline. Important caveat: this is a between-species comparison
(across different loci in different species), not a same-species paired
design - none of the 6 regular-alternation species happened to also have
an identical-copy locus that passed into the control group (4 have no
other core-length group at all; the other 2's identical-copy group
didn't produce a viable independent candidate from tidk).

## Caveats

- **`GAP_THRESHOLD` (50bp) for the telomere-proper boundary is a tunable
  parameter, not validated against ground truth.** For a species with a
  long, genuinely degenerate telomere, a real within-telomere transition
  between variant blocks could in principle exceed 50bp by chance,
  wrongly truncating what's still bona fide telomere. Not tested for
  sensitivity to this choice.
- **Structural verdicts rest on a single genome-wide "best" chromosome-
  end for the formal runs-test z-score** (whichever has the most
  minor-variant hits) - not pooled across all chromosome-ends
  statistically. `max_minor_block`, by contrast, IS a genuine genome-wide
  maximum (checked across every chromosome, both ends, within
  telomere-proper).
- **We tile with the full extracted core (13-15bp for most sets) rather
  than the minimal biological repeat unit (~7bp)**, which is stricter
  than necessary: a single base error anywhere within the core's window
  breaks the match for the entire window, even though it would only
  disqualify one ~7bp period at native resolution. This likely makes our
  block/structure counts throughout this note conservative
  (undercounting true structure), not inflated - not corrected for here.
- **Low-count cases have no statistical power.** Several sets had very
  few minor-variant occurrences even within telomere-proper; their
  "structure" (ISOLATED vs BLOCKY) is only a description of what little
  data exists, not a confirmed statistical conclusion.
- **`Solanum_nigrum`'s core-length-7 set (`AATAAAT`/`ACCTGAA`) involves a
  low-complexity, A/T-rich motif** flagged earlier in this investigation
  as prone to coincidental matches - and indeed, once restricted to
  telomere-proper, it shows *zero* confirmed minor-variant presence (the
  earlier full-window "30.5%, blocky" signal was entirely subtelomeric
  noise). A useful validation of the boundary-restriction fix.
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
- **The identical-copy control group's second-candidate selection is
  inherently noisier than the differing sets' TR-gene-verified variants**
  - candidates come from an unconstrained same-length + small-edit-
  distance search over each species' tidk output, not from an
  independently-validated biological paralog. The decisive filter is
  real telomere-proper presence in the walk (17/62 candidates failed
  this), but a genuinely unrelated repeat family that happens to have
  some real terminal presence for unrelated reasons (e.g. an rDNA-linked
  or satellite repeat near some chromosome ends) could still pass that
  filter without being a true "derivative monomer" of the canonical
  repeat - not fully ruled out for every one of the 62 control sets.

## Reproducing

```bash
# Full (not top-2-only) variant breakdown + rotation-corrected diff positions:
python3 src/tr_core_template_variant_positions.py outputs/tr_repeat_correlation/all_species.tsv \
  outputs/tr_repeat_correlation/core_template_variant_positions.tsv

# Exhaustive telomere-wide dominance + block/runs-test structure over the 44
# differing sets (every main chromosome, both ends; reports both full-window
# and telomere-proper-restricted columns):
python3 src/tr_core_template_array_structure.py outputs/tr_repeat_correlation/core_template_variant_positions.tsv \
  outputs/tr_repeat_correlation/positional 5000 outputs/tr_repeat_correlation/telomere_array_structure.tsv

# Categorize into frequency bands + structure verdicts (z-test-first classification):
python3 src/tr_core_template_array_structure_summary.py outputs/tr_repeat_correlation/telomere_array_structure.tsv \
  outputs/tr_repeat_correlation/telomere_array_structure_summary.tsv

# Control group: find an independent (TR-blind) second candidate for every
# fully-identical-TR-copy species, from that species' own tidk output (also
# writes the variant-positions adapter the next step consumes):
python3 src/find_identical_set_second_candidates.py outputs/tr_repeat_correlation/core_template_sets.tsv \
  outputs/tidk outputs/tr_repeat_correlation/identical_set_second_candidates.tsv \
  outputs/tr_repeat_correlation/identical_set_variant_positions.tsv

# Same exhaustive walk + telomere-proper restriction, applied to the control group:
python3 src/tr_identical_set_array_structure.py outputs/tr_repeat_correlation/identical_set_variant_positions.tsv \
  5000 outputs/tr_repeat_correlation/identical_set_array_structure.tsv

python3 src/tr_core_template_array_structure_summary.py outputs/tr_repeat_correlation/identical_set_array_structure.tsv \
  outputs/tr_repeat_correlation/identical_set_array_structure_summary.tsv

# One-off illustrative render of a specific chromosome-end:
python3 src/walk_terminal_repeat_array.py Centaurium_intermedium SUPER_1_HAP1 3prime CCTAAACC AACCCTAA --window 300
```
