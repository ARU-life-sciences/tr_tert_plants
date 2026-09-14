# Paper planning: headline results, gaps, and file mapping

Snapshot as of 2026-09-13, updated 2026-09-14 with the per-species array-
structure analysis (see #6 below) and again same day after reviewing
directly relevant prior literature (see "Related work" below). Written
to answer: what are the headline results, what's still needed before
this is publication-ready, and which files/notes back each result.

## Related work (read 2026-09-14 - changes framing, not headline validity)

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

**Net effect on framing**: the core mechanism (TR mutation -> telomere
change) and the core pattern (paralog divergence -> heterogeneous
telomere) are now independently demonstrated elsewhere, so the paper's
novelty has to be pitched as generality/scale/synthesis (a
~330-700-species, purely computational test of phenomena previously
shown only in single-genus wet-lab studies, PLUS the first direct link
between specific paralog sequence and specific array structure/block
size - neither prior paper did that combination), not first discovery.

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
   look unless corrected for cyclic rotation. Belyayev et al. 2023
   independently found near-identical block-organized double-monomer
   telomere arrays in Chenopodium (unrelated family), confirmed by
   long-read sequencing and fiber-FISH - but attributed to DNA-level
   recombination, not TR paralogs (untested there). We're the first to
   directly link specific array block structure to specific divergent TR
   paralog sequence - but haven't yet ruled out their mechanism operating
   in our species too (see gap below).

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
- **New, motivated by related work above - two concrete analyses to run
  before this is submittable, both likely cheap with existing data:**
  1. Telomere length vs. paralog-divergence status, dataset-wide (tests
     Kumawat's single-species finding of ~40% longer telomeres in the
     heterogeneous-telomere species, at n=300+) - probably answerable
     from existing `outputs/tidk/` data without new computation.
  2. Whether species with FULLY IDENTICAL TR copies (the 223 identical
     comparable sets, no paralog divergence) ever still show blocky/
     clustered array structure - a direct test of Belyayev's competing
     DNA-recombination-only mechanism against our TR-paralog-linked
     framing. Needs `tr_core_template_array_structure.py`'s machinery run
     against a sample of identical sets instead of only differing ones.

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
