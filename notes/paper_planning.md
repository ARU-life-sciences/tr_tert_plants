# Paper planning: headline results, gaps, and file mapping

Snapshot as of 2026-09-13, updated same day once WEAK_SENSITIVITY cases
were exhaustively resolved (see below). Written to answer: what are the
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
   specific difference is very often physically present in the
   telomere** - confirmed at base-pair resolution across all 44 differing
   sets: 36/44 (82%) show the specific variant genuinely present, exactly
   0/44 are confirmed negatives (the other 8 are the separate LABEL_FLIP
   phenomenon, #6 below). This required exhaustively checking every main
   chromosome's both ends, not a sample - two earlier, smaller-sample
   passes both under-called real positives as absent.
5. **The Template region can be predicted from the gene alone, without
   reference to the telomere** - a boundary-anchored rule (start 2bp
   after the conservation-derived G-rich boundary, take 12bp) recovers
   the true core at mean IoU 0.836 across 1073 loci.

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
- All 8 sensitivity-limited (`WEAK_SENSITIVITY`) cases are now resolved
  (exhaustive check, every main chromosome, both ends): all 8 confirmed
  present. 36/44 is a final number for that bucket, not a floor.
- The 8 `LABEL_FLIP` cases (TR-copy-count "dominant" variant is
  telomerically silent) are unexplained - worth investigating before
  publishing as more than an observation.
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

## Files backing each result

| result | primary file(s) | note documenting it |
|---|---|---|
| Viscum telomerase absence | `outputs/tert_tblastn/Viscum_album.tbl`, `outputs/tr_tbls/Viscum_album*.tbl` (both empty/no confident hit) | *(reported in conversation; not yet its own note - worth writing one)* |
| Carlina non-canonical repeat | `outputs/tr_template_loci/Carlina_vulgaris.tsv`, `outputs/tidk/Carlina_vulgaris.tidk.tsv` | `notes/carlina_vulgaris_telomere_repeat.md` |
| Within-species variation (38% mixed) | `outputs/tr_template_loci/*.tsv` (raw), aggregate via `src/extract_tr_template_all_loci.py` | `notes/within_species_tr_template_variation.md` |
| Fine-scale correlation, 8 species | `outputs/tidk/*`, ad hoc `tidk search` positional scans | `notes/tr_repeat_correlation.md` |
| Corrected core-template stats + 44-set confirmation | `outputs/tr_repeat_correlation/all_species.tsv`, `core_template_sets.tsv`, `final_repeat_confirmations.tsv` | `notes/tr_core_template_variation_and_telomeres.md` |
| Telomere-independent prediction rule | `outputs/tr_repeat_correlation/all_species.tsv` (validation data), `src/predict_tr_template.py` | `notes/tr_template_boundary_prediction.md` |
| Organelle/cobiont data-QC fix (methods caveat, not a finding) | `inputs/dtol_plant_paths.txt`, `src/update_assembly_list.bash` | `notes/organelle_and_cobiont_assembly_pitfall.md` |
