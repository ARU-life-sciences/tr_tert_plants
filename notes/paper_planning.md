# Paper planning: headline results, gaps, and file mapping

Snapshot as of 2026-09-13, updated 2026-09-14 with the per-species array-
structure analysis (see #6 below). Written to answer: what are the
headline results, what's still needed before this is publication-ready,
and which files/notes back each result.

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
   statistic).
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
   it clusters** - 41/44 (93%) of differing sets show contiguous,
   block-like structure (2-299 units), not scattered single-base
   substitutions; 5/44 show a striking regular short-period alternation
   between two variants (a candidate genuine higher-order-repeat, HOR,
   signature - e.g. `Centaurium_intermedium` alternates almost perfectly
   for 18+ consecutive units). Also: 25/44 sets have a telomere-dominant
   sequence that differs from the TR-copy-count dominant (broader than
   the original 8 LABEL_FLIP cases - only found once tiling used every
   distinct variant, not just the top 2). And a real methodological catch
   along the way: 6/55 "differing" variants across the dataset are pure
   phase reads of the identical periodic repeat (0 true substitutions),
   not point-mutant paralogs at all - naive position-wise comparison of
   these short, near-periodic motifs overstates how different two cores
   look unless corrected for cyclic rotation.

## What's needed before this is paper-ready

**Cross-cutting, applies to all findings:**
- No literature search has been done at all. Needs a proper comparison
  against Fajkus et al. 2019 (the paper this whole domain-map is built
  on) and a check for prior reports of Viscum/mistletoe telomerase loss
  or Asteraceae telomere variation.
- Everything is single-haplotype, single-assembly. No biological
  replication.
- No wet-lab validation anywhere (PCR, TRAP assay, Southern/FISH) -
  everything is computational.
- Decide manuscript scope: this is plausibly 2-3 papers (Viscum loss;
  Carlina repeat-switch; the dataset-wide fine-scale correlation as a
  resource/methods paper), not one.

**Specific to #1 (Viscum):**
- Check assembly completeness/QC in the relevant regions - need to rule
  out "assembly gap" as an alternative to "genuinely absent."
- No phylogenetic context (nearby parasitic/non-parasitic relatives) has
  been added yet.

**Specific to #2 (Carlina):**
- Which locus (if any) is the *functional* copy is unresolved - 6 loci
  show the variant, 1 weak locus looks canonical.
- Second haplotype not checked.
- Direct/reverse-complement orientation ambiguity (see below) directly
  affects whether the templating claim is stated correctly.

**Specific to #3/#4 (dataset-wide):**
- **Haplotype-phasing was never systematically screened** -
  `Empetrum_nigrum`'s `SUPER_N_HAP2/3/4` naming showed this is a real
  risk, not hypothetical; some "independent loci" could be allelic
  copies of one gene, inflating copy-number claims. This needs a proper
  filter before any copy-number statistic goes in a paper.
- All 16 sensitivity-limited cases (8 `WEAK_SENSITIVITY` + 8
  `LABEL_FLIP`) are now resolved (exhaustive check, every main
  chromosome, both ends): all 16 confirmed present. 44/44 is a final
  number, not a floor.
- **Orientation ambiguity is unresolved and matters**: 85% of matches
  are in the raw `direct` sense rather than the textbook
  reverse-complement templating relationship. A reviewer will ask about
  this; it needs an answer or an explicit acknowledged limitation.
- Multiple-testing correction across the ~330-species, ~1450-locus
  screen is informal (per-locus shuffle test only) - should be
  formalized.

**Specific to #5 (prediction rule):**
- Validated only against its own discovery dataset - no held-out test
  set. Needs validation on species/loci not used to derive the rule.

**Specific to #6 (array structure):**
- Structural verdicts (blocky vs regular-alternation) rest on a single
  best-powered chromosome-end per set for the formal runs-test z-score,
  not a statistic pooled across all chromosome-ends - needs a proper
  pooled/multi-chromosome test before a p-value goes in a paper.
- The regular short-period alternation (candidate HOR) in 5 sets hasn't
  been checked against any independent expectation (e.g. does the
  alternating unit correspond to anything at the sequence/structural
  level, or could it be an artifact of how two similar-length candidate
  strings compete for the same genomic slots) - needs closer inspection
  before calling it a genuine HOR.
- 4 of the "differing" sets turned out to be pure phase-rotations of one
  repeat, not real sequence differences - their inclusion in earlier
  copy-number/diversity statistics (headline #3) hasn't been corrected
  for; worth rechecking whether removing them changes those percentages
  materially (likely marginal at n=44, but should be checked, not
  assumed).

## Files backing each result

| result | primary file(s) | note documenting it |
|---|---|---|
| Viscum telomerase absence | `outputs/tert_tblastn/Viscum_album.tbl`, `outputs/tr_tbls/Viscum_album*.tbl` (both empty/no confident hit) | *(reported in conversation; not yet its own note - worth writing one)* |
| Carlina non-canonical repeat | `outputs/tr_template_loci/Carlina_vulgaris.tsv`, `outputs/tidk/Carlina_vulgaris.tidk.tsv` | `notes/carlina_vulgaris_telomere_repeat.md` |
| Within-species variation (38% mixed) | `outputs/tr_template_loci/*.tsv` (raw), aggregate via `src/extract_tr_template_all_loci.py` | `notes/within_species_tr_template_variation.md` |
| Fine-scale correlation, 8 species | `outputs/tidk/*`, ad hoc `tidk search` positional scans | `notes/tr_repeat_correlation.md` |
| Corrected core-template stats + 44-set confirmation | `outputs/tr_repeat_correlation/all_species.tsv`, `core_template_sets.tsv`, `final_repeat_confirmations.tsv` | `notes/tr_core_template_variation_and_telomeres.md` |
| Telomere-independent prediction rule | `outputs/tr_repeat_correlation/all_species.tsv` (validation data), `src/predict_tr_template.py` | `notes/tr_template_boundary_prediction.md` |
| Per-species array structure (clustering, true dominance, HOR candidates) | `outputs/tr_repeat_correlation/telomere_array_structure.tsv`, `telomere_array_structure_summary.tsv`, `core_template_variant_positions.tsv` | `notes/tr_core_template_array_structure.md` |
| Organelle/cobiont data-QC fix (methods caveat, not a finding) | `inputs/dtol_plant_paths.txt`, `src/update_assembly_list.bash` | `notes/organelle_and_cobiont_assembly_pitfall.md` |
