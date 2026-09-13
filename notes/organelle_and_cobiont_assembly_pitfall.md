# Pitfall: "latest assembly" can silently pick a plastid or cobiont genome

## Summary

Two species in the DToL plant assembly tree have more than one curated
assembly directory - a real nuclear genome, plus one or more
plastid/cobiont assemblies from the same sequencing run - and in both
cases the non-nuclear assembly happened to be curated *more recently*.
`any_latest_assemblies.bash` (the tool `src/update_assembly_list.bash`
wraps) picks whichever directory has the newest mtime with no awareness of
assembly type, so it silently picked the wrong genome for both:

- **Galanthus_nivalis**: only a plastid assembly exists at all right now
  (`lsGalNiva1.plastid.1`, curated Jul 2026) - the nuclear genome
  (`lsGalNiva1.1`) is still mid-pipeline (contamination/QC filtering only,
  no finished primary fasta yet). Removed from `inputs/dtol_plant_paths.txt`
  entirely until DToL publishes a real nuclear assembly.
- **Lemna_minuta**: has a real nuclear genome (`laLemMinu1.1`, curated Nov
  2024, 89-90MB), but also ~40 separately-curated *microbiome cobiont*
  assemblies (`laLemMinu1.Leptothrix_sp_1.1`, `laLemMinu1.Acinetobacter_sp_1.1`,
  etc. - duckweed's associated bacteria, sequenced and assembled alongside
  the host), all curated later (Apr 2025). The picker grabbed
  `Leptothrix_sp_1` (~1MB, a single bacterium) instead of the actual plant
  genome. Fixed to point at `laLemMinu1.1` and every stage re-run.

Both were, until caught, silently reported as legitimate biological
findings ("zero TERT hits") when checked casually - they *looked* like
Viscum_album's genuine result (same shape: clean job completion, empty
output), which is exactly what made them dangerous. They were only caught
because a multi-locus variation check (unrelated goal) turned up a
suspiciously tiny genome size for Galanthus, which prompted checking the
path, which led to checking the whole dataset for the same pattern.

**Viscum_album's result is unaffected** - it has only one curated
directory, correctly named, and its assembly is genome-scale (~40GB peak
memory during the tblastn search, matching a real nuclear genome, not a
~1-150KB organelle/single-bacterium assembly).

## Fix

`src/update_assembly_list.bash` now screens every path from the "latest
assembly" tool before it's allowed into `inputs/dtol_plant_paths.txt`:

1. **Naming check**: rejects any path whose curated-directory name doesn't
   match the plain `PREFIX<digits>(.hap<digits>)?.<digits>` pattern (e.g.
   `daCarVulg1.2`, `drDryOcto1.hap1.2`) - anything with an extra
   organism-name segment (`.plastid.`, `.Leptothrix_sp_1.`, `.metagenome.`,
   etc.) is flagged.
2. **Size floor**: as a backstop for cases that don't announce it in the
   directory name, rejects any assembly file under 5MB compressed (`MIN_BYTES`,
   overridable) - comfortably below any real nuclear plant genome, well
   above a compressed plastid or single-bacterium assembly.

Flagged paths are excluded from the fresh list (not silently added, not
silently updated) and reported separately for manual investigation - see
`update_assembly_list.bash`'s `SUSPECT` output.

Checked the naming pattern across the other ~700 species already in the
list: nothing else currently matches it, so this looks like an isolated
(if serious) issue rather than a systemic one across the dataset - but
it's now a permanent, automatic check on every future refresh rather than
something that has to be noticed by chance again.

## Caveat

This screen is heuristic, not exhaustive. A cobiont assembly that happens
to be large (e.g. a fungal endophyte, which can be tens of MB) and
happens to follow the plain naming pattern (unlikely given DToL's
observed convention of naming cobiont directories after the organism, but
not guaranteed) could still slip through. If a species' TERT/TR result
looks like a clean, total absence, it's worth a quick sanity check of the
genome path and assembly size before treating it as a finding - the same
check that caught these two.
