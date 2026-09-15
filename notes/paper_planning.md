# Paper planning: headline results, gaps, and file mapping

Snapshot as of 2026-09-13, updated 2026-09-14 with the per-species array-
structure analysis (see #6 below), again after reviewing directly
relevant prior literature (see "Related work" below), and again after
adding a controlled TR-identical-copy comparison group plus fixing two
real methodology issues in the array-structure analysis (subtelomeric
noise, a verdict-classification bug - see #6). Updated again 2026-09-16
after resolving two more cross-cutting gaps: the haplotype-phasing
redundancy screen (real but narrow - 8/331 species affected, only 1 of
the 44 headline differing sets) and the orientation-ambiguity question
(resolved as a strand-convention artifact, not a biological anomaly -
see "Specific to #3/#4" below), plus the multiple-testing correction,
pure-rotation correction, and `GAP_THRESHOLD` sensitivity check (all
DONE, see #6/#3/#4 below). Updated again 2026-09-15: manuscript scope
decided (one paper, not 2-3), and a major new related-work paper found
(Závodník et al. - see below, materially revises the novelty framing).
Written to answer: what are the headline results, what's still needed
before this is publication-ready, and which files/notes back each
result.

## Related work (read 2026-09-14, 2026-09-15 - changes framing, not headline validity)

- **Fajkus et al. 2019, NAR** ("Telomerase RNAs in land plants") - the
  paper our own TR HMM is seeded from (confirmed: `inputs/TR.fasta`
  matches their identified orthologs, including the corrected `AtTR`,
  not the disproven `TER1`). Directly relevant: (i) template region
  defined as repeat-unit-length + >=1nt, matching why our own core
  lengths cluster there; (ii) **single-nucleotide TR template mutations
  in vitro produce exactly the corresponding telomere repeat change** -
  direct wet-lab proof of the mechanism our whole analysis assumes; (iii)
  already predicted mixed-telomere-motif TRs in some Asparagales species,
  i.e. our headline #3/#4 phenomenon was suspected there, just not
  confirmed base-pair-resolved or at scale.
- **Kumawat et al. 2025, PLOS Genetics** (Mimulus) - a near-direct
  single-genus precedent: TR duplicates, paralogs diverge in template
  sequence, retained-and-expressed paralogs produce a genuinely
  heterogeneous telomere (confirmed via TRF, FISH, and Sanger-sequenced
  TRAP clones - *M. lewisii*, 9 vs 26 reads for its two repeat types);
  a paralog can be expressed but barely used (*M. cardinalis*, ~9.7%
  minor frequency - in our own "RARE" band range); duplication -> divergence
  -> retain-if-functional/pseudogenize-if-not is their proposed model.
  They also found the heterogeneous-telomere species has ~40% longer
  telomeres than its single-repeat sister species.
- **Belyayev et al. 2023, BMC Genomics** (Chenopodium) - independently
  describes "block-organized double-monomer terminal telomeric arrays"
  (blocks of canonical repeat interchanging with a derivative monomer),
  confirmed via long-read sequencing, fiber-FISH, and a synthetic
  junction-spanning probe - close to our own blocky/regular-alternation
  structural classes (headline #6), in an unrelated plant family. Their
  proposed mechanism is DNA-level recombination between G-/C-rich
  strands, NOT TR-paralog-driven (Chenopodium wasn't checked for TR
  duplication) - a genuine competing hypothesis we haven't ruled out.
- **Závodník et al.** (read 2026-09-15; exact journal/year not yet
  confirmed - shares methodology/covariance model with Fajkus et al.
  2021, likely same group) - the single most directly relevant paper
  found so far, more so than Kumawat or Belyayev. Screened ~986
  Tracheophyta genomes with Infernal/covariance models (comparable raw
  scale to our own ~330-700 species) and built a large table of species
  with multiple TR paralogs whose templates predict different telomere
  motifs, validating a curated subset (~20-40 species) with FISH,
  BAL31-TRF, and TRAP+sequencing. Their proposed model (Fig. 1): TR
  duplication -> relaxed selection on the redundant copy -> mutation
  accumulates -> a transitional "mixed telomere" phase (both ancestral
  and mutant motifs coexisting) -> resolved by either selection
  (pseudogenize/repair the mutant, revert to ancestral) or genetic drift
  (mutant motif replaces ancestral, completing a transition) - a more
  complete, directly citable theoretical frame than Kumawat's simpler
  retain-if-functional/pseudogenize-if-not dichotomy; worth adopting as
  our interpretive lens. Most relevant to headline #6: FISH typing of
  mixed-telomere species into three spatial patterns (high-colocalization
  "mixed", low-colocalization/partial segregation, single-motif), and -
  independently, via PacBio long reads, a completely different sequencing
  technology from ours - *Fagopyrum tataricum*'s two telomere motifs form
  "intermingled" arrays with **no homogeneous stretches of either motif**,
  essentially the same phenomenon as our regular-alternation/candidate-HOR
  finding. Also links TR *transcript abundance* (RT-qPCR), not just gene
  copy number, to which motif dominates the real telomere in Fagopyrum -
  a more mechanistically direct version of our "TR-copy-count doesn't
  predict telomere dominance" finding.

**Net effect on framing - revised again after Závodník et al.**: raw
scale is no longer a strong novelty claim on its own (Závodník et al.
already screened a comparable number of genomes). What remains
genuinely distinctive: (1) a fully computational, exhaustive,
base-pair-resolved confirmation across *every* chromosome end for
*every* species with paralog divergence in our dataset, vs their
curated wet-lab validation of ~20-40 selected species; (2) critically,
**a formally statistically-tested comparison against a matched
TR-identical control group** (our Fisher's-exact 6/44-vs-0/62 regular-
alternation result) - nothing like this exists in any of the four
related papers, whose structural characterizations are qualitative/
descriptive, case-by-case, not a dataset-wide hypothesis test against a
proper null. Lead the paper with the controlled-comparison result, not
with sample size.

## Headline results

1. **`Viscum_album` (mistletoe) shows convergent absence of both
   telomerase components** — no confident TERT hit (tblastn, including
   one exhaustive 12.4-day search) and no confident TR hit (nhmmer) —
   candidate genuine telomerase loss in a parasitic plant.
2. **`Carlina_vulgaris` carries a candidate independently-evolved,
   non-canonical telomere repeat** (`AAACTG`-family, not `TTTAGGG`) —
   supported by 6/7 TR gene paralogs *and* independent genome-wide
   positional confirmation (real chromosome-terminal enrichment for the
   variant, none for canonical).
3. **TR gene paralog diversity within a single genome is common, not
   rare**: 38% of multi-locus species carry both canonical- and
   divergent-Template copies simultaneously; at the properly-resolved
   templating core specifically, 71-84% of copies are identical even at
   high copy number (correcting an initial inflated 97% pairwise
   statistic). Kumawat et al. 2025 showed the same paralog-duplication ->
   divergence pattern in one genus (Mimulus, wet-lab validated); our
   contribution is showing it's general across land plants at
   dataset-wide scale, not a Mimulus peculiarity.
4. **Where TR gene copies genuinely differ at the templating core, that
   specific difference is essentially always physically present in the
   telomere** - confirmed at base-pair resolution across all 44 differing
   sets: 44/44 (100%) show the specific variant genuinely present, 0/44
   are confirmed negatives. This required exhaustively checking every
   main chromosome's both ends, not a sample - two earlier, smaller-sample
   passes both under-called real positives as absent, and this also
   overturned the original "LABEL_FLIP" framing (8 sets initially read as
   the TR-copy-majority variant being telomerically silent) - it isn't
   silent in any of the 8, just present at very low relative frequency
   (0.1-3.7% in 7/8 cases; the 8th, `Solanum_nigrum` at 30.5%, involves a
   low-complexity motif and deserves extra scrutiny). TR gene paralog
   copy number predicts telomeric *abundance* of a variant, not its
   presence/absence.
5. **The Template region can be predicted from the gene alone, without
   reference to the telomere** - a boundary-anchored rule (start 2bp
   after the conservation-derived G-rich boundary, take 12bp) recovers
   the true core at mean IoU 0.836 across 1073 loci.
6. **The minor variant's arrangement within the telomere is not random -
   it clusters - and a controlled comparison shows regular short-period
   alternation (candidate HOR) specifically requires TR paralog
   divergence.** After two real methodology fixes (a fixed 5000bp window
   can include subtelomeric/interstitial sequence past the true telomere,
   now restricted to "telomere-proper"; and a verdict-classification bug
   that silently absorbed regular-alternation cases into a generic
   "clustered" label, now fixed by checking the formal runs-test
   statistic first): 35/41 (85%) of differing sets with a confirmed minor
   variant show non-random structure. The key new result: run the
   identical method on 62 sets from species with **fully identical** TR
   copies (no paralog divergence, using an independently tidk-discovered
   candidate blind to TR content) as a control group. Ordinary
   blocky/clustered structure happens at a real baseline rate regardless
   of TR status (30/62, consistent with Belyayev et al. 2023's proposed
   DNA-recombination mechanism operating generally) - but **regular
   short-period alternation was found in 6/44 differing sets and 0/62
   controls (Fisher's exact p=0.004)**, and any confirmed minor variant
   at all is significantly more common in differing sets (93.2% vs 72.6%,
   p=0.011) - both robust to a 2x change in the telomere-proper boundary
   parameter (checked 2026-09-16; a third, weaker sub-comparison was not
   robust and should be treated as suggestive only, not reported with the
   same confidence). This is the first direct evidence tying the
   fine-grained alternating pattern specifically to TR paralog
   divergence, not just generic recombination noise - something none of
   Kumawat et al., Belyayev et al., or Závodník et al. could test, since
   none had both a TR-divergent and TR-identical group to compare
   (Závodník et al. independently found the same qualitative pattern -
   PacBio-confirmed "intermingled," non-homogeneous dual-motif arrays in
   *Fagopyrum tataricum* - via a different sequencing technology, which
   is real, useful corroboration of the phenomenon itself, but their
   characterization stays descriptive/case-by-case, not a statistical
   test against a null).
   Also: 25/44 sets have a telomere-dominant sequence that differs from
   the TR-copy-count dominant. And a real methodological catch along the
   way: 6/55 "differing" variants across the dataset are pure phase reads
   of the identical periodic repeat (0 true substitutions), not
   point-mutant paralogs at all - naive position-wise comparison of these
   short, near-periodic motifs overstates how different two cores look
   unless corrected for cyclic rotation.

## What's needed before this is paper-ready

**Cross-cutting, applies to all findings:**
- Literature search started 2026-09-14 (see "Related work" above:
  Fajkus 2019, Kumawat 2025, Belyayev 2023) - still needed: Viscum/
  mistletoe telomerase loss precedent, Asteraceae telomere variation
  beyond the two Asteraceae sources already found (an old unpublished
  Garnatje grant report specifically flags *Arctium* clustering - matches
  our own `Arctium_lappa`/`minus` results - and Mlinarec et al. 2019 on
  *Tanacetum* subtelomeric repeats), and the Závodník et al. paper
  Kumawat cites for TR duplication across multiple plant families.
- Everything is single-haplotype, single-assembly. No biological
  replication.
- No wet-lab validation anywhere (PCR, TRAP assay, Southern/FISH) -
  everything is computational.
- **DECIDED (2026-09-15): one paper**, not a 2-3 split. Viscum is
  included as "what we have" (ongoing work elsewhere, not pursued
  further here - see #1's gaps below) rather than built out into its
  own paper; Carlina and the dataset-wide findings are folded into the
  same manuscript.

**Specific to #1 (Viscum):**
- **Not being pursued further here** (per the one-paper scope decision
  above) - Viscum work is ongoing elsewhere; this note reports what we
  already have (no confident TR or TERT hit) as-is, not as a fully
  worked-up standalone finding. The gaps below are real but out of scope
  for now:
  - Check assembly completeness/QC in the relevant regions - need to
    rule out "assembly gap" as an alternative to "genuinely absent."
  - No phylogenetic context (nearby parasitic/non-parasitic relatives)
    has been added yet.

**Specific to #2 (Carlina):**
- Which locus (if any) is the *functional* copy is unresolved - 6 loci
  show the variant, 1 weak locus looks canonical.
- Second haplotype not checked.
- Direct/reverse-complement orientation question is now resolved (see
  below) - doesn't block the templating claim.

**Specific to #3/#4 (dataset-wide):**
- **DONE (2026-09-16): haplotype-phasing screen.** `Empetrum_nigrum`'s
  `SUPER_N_HAP2/3/4` naming flagged the risk; a full sweep (covering both
  DToL haplotype-naming conventions, `SUPER_N_HAPx` and
  `HAPx_SCAFFOLD_N`/`HAPx_SUPER_N`) found 68/331 species have >=1 TR
  locus on a haplotype-labelled scaffold, but only **8/331** show the
  genuinely inflating pattern (same true chromosome, TR loci on >1
  haplotype number): `Buddleja_davidii`, `Empetrum_nigrum`,
  `Galium_boreale`, `Hesperis_matronalis`, `Hypericum_perforatum`,
  `Lythrum_salicaria`, `Salix_cinerea`, `Solidago_canadensis`. In every
  case the redundant haplotype copies carry identical sequence (allelic
  redundancy of one real locus, not fabricated diversity) - inflates raw
  copy-number counts for these 8 but doesn't create spurious "differing"
  classifications. Only 1 of the 44 headline differing sets is affected
  (`Lythrum_salicaria`, 11bp - true locus count is 3, not 6; see the
  corrected table in `tr_core_template_variation_and_telomeres.md`).
  Script: `src/check_haplotype_redundancy.py`, output:
  `outputs/tr_repeat_correlation/haplotype_redundancy_flagged.tsv`.
- All 16 sensitivity-limited cases (8 `WEAK_SENSITIVITY` + 8
  `LABEL_FLIP`) are now resolved (exhaustive check, every main
  chromosome, both ends): all 16 confirmed present. 44/44 is a final
  number, not a floor.
- **DONE (2026-09-16): orientation ambiguity resolved.** 85% of matches
  are in the raw `direct` sense rather than the textbook
  reverse-complement templating relationship - checked, and it's a
  strand-convention artifact, not a biological anomaly: 98.8% of
  `direct` matches are C-rich Template sequences matching a C-rich
  `tidk`-reported repeat. Both are independently correctly C-rich (our
  Template extraction matches Fajkus et al. 2019's own validated C-rich
  `AtTR` template; `tidk`'s reported repeat is C-rich most plausibly
  because `AAACCCT` is the global lexicographic minimum among all 14
  rotations of the canonical repeat family, a standard convention for
  representing a strand-ambiguous repeat). No confound with which
  genomic strand the TR gene itself sits on (~50/50 regardless of match
  orientation, checked across all 1450 loci). See
  `tr_repeat_correlation.md`'s caveats for the full check.
- **DONE (2026-09-16): multiple-testing correction.** BH-FDR across the
  1293 tested loci: 1083/1293 (83.8%) significant at q<0.05, 990/1293
  (76.6%) at q<0.01 - comfortably supports the significance claim (more
  loci pass formal FDR than the naive "exactly 0/200 shuffles" count).
  One real remaining limitation: q<0.001 is unreachable at only 200
  shuffles/locus (floor p~0.005) - would need more shuffles for a
  stricter published threshold. Script:
  `src/multiple_testing_correction.py`.

**Specific to #5 (prediction rule):**
- Validated only against its own discovery dataset - no held-out test
  set. Needs validation on species/loci not used to derive the rule.

**Specific to #6 (array structure):**
- **DONE (2026-09-14): identical-TR-copy control group.** 62 sets from
  60 species with fully identical TR copies, walked with the same
  exhaustive method using an independently tidk-discovered candidate.
  Regular alternation: 0/62 controls vs 6/44 differing sets (p=0.004).
  This is a between-species comparison, not a same-species paired design
  - none of the 6 regular-alternation species had a matched
  identical-copy locus that passed into the control group. A same-species
  paired test (if a species with BOTH a differing and a passing identical
  set turns up in future data) would be stronger evidence still.
- **Telomere length vs. paralog-divergence: not attempted.** Assembled
  "telomere length" from draft genome assemblies conflates real biology
  with assembler/sequencing-technology choices about how much of a
  repetitive terminal array to resolve; we have no assembly-quality
  metadata (N50, technology, curation batch) to control for that
  confound, and no wet-lab length validation. Judged not reliably
  answerable with current data - dropped rather than reported with an
  uninterpretable result.
- Structural verdicts (blocky/clustered vs regular-alternation) still
  rest on a single best-powered chromosome-end per set for the formal
  runs-test z-score, not a statistic pooled across all chromosome-ends -
  needs a proper pooled/multi-chromosome test before a p-value goes in a
  paper.
- **DONE (2026-09-16): `GAP_THRESHOLD` sensitivity check.** Reran both
  groups at 100bp instead of the default 50bp. Two of the three headline
  comparisons are robust: regular alternation is **exactly** 6/44 vs
  0/62 at both thresholds (p=0.004 both times); "any confirmed minor
  variant" stays significant (if anything slightly stronger: p=0.004 at
  100bp vs p=0.011 at 50bp). The third ("non-random structure, given
  present") is **not** robust - already-borderline p=0.049 at 50bp
  becomes p=0.42 (not significant) at 100bp. Report the first two with
  confidence; drop or clearly caveat the third if it goes in a paper.
- The regular short-period alternation (candidate HOR) in 6 sets still
  hasn't been checked against any independent expectation at the
  sequence/structural level - needs closer inspection before calling it
  a genuine HOR, though its exclusive association with TR-paralog
  divergence (vs. 0/62 controls) is itself real evidence against a pure
  mapping-artifact explanation.
- **DONE (2026-09-16): pure-rotation correction.** Reclassifying the 4
  pure-phase-rotation sets as identical (not genuinely divergent) moves
  headline #3's conservation stats up modestly: 83.5%->85.0% overall,
  71.1%->73.7% for the >=5-loci subset. Also drops the differing-set
  count from 44/33 species to 40/30 species (`Platanus_x_hispanica`,
  `Acaena_novae_zelandiae`, `Acaena_ovalifolia` are dropped entirely -
  each had every one of its differing sets turn out to be pure rotation).
  Confirmed marginal, as expected, not headline-changing. Script:
  `src/pure_rotation_corrected_stats.py`.
- The identical-copy control group's second-candidate selection is
  inherently noisier than the differing sets' TR-gene-verified variants -
  see the note's caveats section.
- We tile with the full extracted core (13-15bp) rather than the minimal
  biological repeat unit (~7bp), which is a stricter match criterion than
  necessary and likely makes block/structure counts throughout this
  finding conservative (undercounts), not inflated - not corrected for.

## Files backing each result

| result | primary file(s) | note documenting it |
|---|---|---|
| Viscum telomerase absence | `outputs/tert_tblastn/Viscum_album.tbl`, `outputs/tr_tbls/Viscum_album*.tbl` (both empty/no confident hit) | *(reported in conversation; not yet its own note - worth writing one)* |
| Carlina non-canonical repeat | `outputs/tr_template_loci/Carlina_vulgaris.tsv`, `outputs/tidk/Carlina_vulgaris.tidk.tsv` | `notes/carlina_vulgaris_telomere_repeat.md` |
| Within-species variation (38% mixed) | `outputs/tr_template_loci/*.tsv` (raw), aggregate via `src/extract_tr_template_all_loci.py` | `notes/within_species_tr_template_variation.md` |
| Fine-scale correlation, 8 species | `outputs/tidk/*`, ad hoc `tidk search` positional scans | `notes/tr_repeat_correlation.md` |
| Corrected core-template stats + 44-set confirmation | `outputs/tr_repeat_correlation/all_species.tsv`, `core_template_sets.tsv`, `final_repeat_confirmations.tsv` | `notes/tr_core_template_variation_and_telomeres.md` |
| Telomere-independent prediction rule | `outputs/tr_repeat_correlation/all_species.tsv` (validation data), `src/predict_tr_template.py` | `notes/tr_template_boundary_prediction.md` |
| Per-species array structure (clustering, true dominance, HOR candidates) + TR-identical control group | `outputs/tr_repeat_correlation/telomere_array_structure.tsv`, `telomere_array_structure_summary.tsv`, `core_template_variant_positions.tsv`, `identical_set_second_candidates.tsv`, `identical_set_array_structure.tsv`, `identical_set_array_structure_summary.tsv` | `notes/tr_core_template_array_structure.md` |
| Organelle/cobiont data-QC fix (methods caveat, not a finding) | `inputs/dtol_plant_paths.txt`, `src/update_assembly_list.bash` | `notes/organelle_and_cobiont_assembly_pitfall.md` |
