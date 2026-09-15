# Epilobium hirsutum: candidate TERT loss despite an intact, functional TR

**Status: candidate finding, comparative-genomic evidence only.** No
wet-lab confirmation; not checked against the literature for this
species/genus.

## Summary

*Epilobium hirsutum* (great willowherb, Onagraceae) shows no confident
TERT (telomerase reverse transcriptase) hit by any of 3 HMM profiles or
by tblastn - the same pattern as `Viscum_album` (headline #1) - but
unlike Viscum, its TR (telomerase RNA) gene is clearly present (3 loci,
all sharing an identical core template, `AAACCCTAAACCT`) and was already
independently confirmed elsewhere in this investigation (the
TR-identical control group built for headline #6) to genuinely match
real telomeric sequence in an exhaustive, base-pair-resolved,
telomere-proper-restricted walk - sparse (3 confirmed occurrences
dataset-wide) but real and non-zero. A targeted synteny check against
its congener *E. ciliatum*
(which has a strong, unambiguous TERT hit) rules out assembly
incompleteness as the explanation: the genomic neighbourhood around
TERT maps cleanly and confidently to a specific, fully-resolved,
mid-scaffold locus in *E. hirsutum* - but the TERT gene body itself has
no detectable sequence similarity there at all, even at permissive
nucleotide-level blastn stringency.

This is a materially different, and arguably more interesting, case
than Viscum: a TR-copy-retaining, apparently-still-functional-templating
species with a specific, localised TERT absence, rather than the
simpler "both components gone" pattern.

## Background

The standard pipeline searches every DToL species for TERT via 3 HMM
profiles (`catalytic_subunit`, `isoformx1`, `tert1` - `run_tert_hmms.bash`)
followed by tblastn confirmation against curated plant TERT references
(`run_tblastn_tert_bsub.bash`). `Epilobium_hirsutum` returns zero hits
from all four searches:

```
outputs/tert_tbls/Epilobium_hirsutum/*.tbl   - all 3 HMM profiles: 0 hits
outputs/tert_tblastn/Epilobium_hirsutum.tbl  - 0 bytes
```

Its TR gene search (`outputs/tr_tbls/Epilobium_hirsutum.tbl`), by
contrast, is a normal, non-empty, multi-locus hit: 3 loci, all sharing
one identical 13bp core template (`AAACCCTAAACCT`) -
`core_template_sets.tsv` row `Epilobium_hirsutum 13 3 1 AAACCCTAAACCT`.
This made it eligible for, and it was included in, the TR-identical
control group built for headline #6
(`outputs/tr_repeat_correlation/identical_set_array_structure.tsv`),
where the same exhaustive, telomere-proper-restricted walk confirmed the
TR core genuinely present in the real telomere (`proper_dom_n=3` - a low
but non-zero, real confirmed count; the candidate second/minor variant
tested alongside it in that analysis was not confirmed present, unrelated
to the point here). So the TR side of telomerase is not in question;
only TERT.

The obvious alternative explanation for "no hit" is assembly
incompleteness (the same caveat already flagged, unresolved, for
Viscum), not genuine absence - draft genome assemblies can have real
gaps, and a single-copy gene sitting in one is indistinguishable from a
genuinely deleted one by a simple hit/no-hit search alone. This needed a
positive, independent check, not just "we didn't find it."

## Evidence: synteny check against a congener with a confident TERT hit

`Epilobium_ciliatum` (same genus) has an unambiguous, strong TERT hit -
consistent multi-exon matches across all 3 HMM profiles and tblastn, at
e-values down to 1e-36, all localised to one region of `SUPER_16`
(~2,602,800-2,609,900). That rules out "our TERT reference library
doesn't cover this lineage" - a close relative is trivially found.

`src/check_gene_synteny_gap.py` extracts a 20kb window around this
region in `E. ciliatum` (`SUPER_16:2,595,000-2,615,000`, generously
padded on both sides of the gene) and blastn's it against `E. hirsutum`'s
full genome (e-value <= 1e-10):

```bash
python3 src/check_gene_synteny_gap.py Epilobium_ciliatum SUPER_16 2595000 2615000 \
  Epilobium_hirsutum outputs/synteny_checks
```

Result:

- **A confident, single-locus anchor exists**: query positions
  3,794-6,296 (i.e. `E. ciliatum` genome coordinates ~2,598,800-2,601,300,
  just upstream of the TERT gene) match `E. hirsutum SUPER_16:3,043,979-
  3,047,576` at 87-93% identity, colinear, best e-value 0.0 - a real,
  single-locus syntenic anchor, not scattered repetitive-element noise
  (repetitive-element hits, by contrast, scatter across a dozen-plus
  different scaffolds at once - seen clearly flanking this anchor on
  both sides in the raw hit table, and correctly not treated as
  informative synteny).
- **The TERT gene body itself has zero hits**: query positions
  6,883-15,353 (~8.5kb, comfortably spanning the entire TERT ORF region
  based on `E. ciliatum`'s own coordinates) returns **no blastn hit
  anywhere in `E. hirsutum`'s genome**, at any identity, even at the
  permissive 1e-10 threshold used for the whole search.
- **Not an assembly gap at this locus**: the anchor sits on
  `SUPER_16`, a 10,471,454bp scaffold, at position ~3.04Mb - 29% of the
  way through, nowhere near either end. A 30kb window centred on the
  anchor (`SUPER_16:3,030,076-3,060,076`) contains **0% N's** - fully
  resolved, real, sequenced DNA, not unresolved/gapped assembly.

So the genomic neighbourhood is present, confidently localised to one
specific place, and completely resolved - but whatever sequence
actually occupies the interval where TERT should be, it isn't
detectably TERT, at a locus where the same search finds TERT easily in
a close relative sharing 87-93% flanking identity.

## Interpretation

The evidence is consistent with a genuine, localised loss (deletion, or
divergence far beyond what the flanking sequence shows) of the TERT
gene body specifically in *E. hirsutum*, with an otherwise intact and
correctly assembled genomic neighbourhood - not an assembly-completeness
artifact. Combined with a TR gene that's present, multi-copy, and
already confirmed to genuinely template the real telomere, this would
be a case of retained/functional TR with an apparently missing TERT -
the reverse-transcriptase partner, not the templating RNA, being the
lost component. Whether this reflects true telomerase loss (with some
telomerase-independent maintenance mechanism substituting, as documented
in other eukaryote lineages - see `tr_repeat_correlation.md`'s
background) or an as-yet-undetected, highly diverged TERT paralog
elsewhere in the genome isn't resolved here.

## Caveats

- **Single-congener comparison only.** One relative, one direction of
  comparison. Doesn't rule out a large-scale genomic rearrangement that
  moved (rather than deleted) TERT to a completely different, unlinked
  part of the `E. hirsutum` genome, which a local synteny window
  wouldn't detect - a whole-genome low-stringency TERT search (rather
  than the standard pipeline's curated-reference tblastn) has not been
  run to check for this.
- **Only one haplotype checked, and the two compared species use
  different ones**: `E. hirsutum`'s assembly is `drEpiHirs1.hap2.1.primary`
  (haplotype 2) while `E. ciliatum`'s is `drEpiCili1.hap1.1.primary`
  (haplotype 1). Shouldn't matter for a species-level presence/absence
  question, but worth confirming `E. hirsutum`'s hap1 (if available)
  agrees before treating this as settled.
- **No wet-lab validation** - this is purely comparative-genomic
  evidence, the same evidentiary tier as the Viscum finding it's being
  compared against, not confirmed by an independent method.
- **The TERT-body "gap" interval (~8.5kb) is defined relative to
  `E. ciliatum`'s own gene coordinates**, not an independently annotated
  `E. hirsutum` gene model (no annotation was available) - the exact
  boundaries of what's missing vs merely diverged aren't pinned down at
  base-pair resolution.
- **No literature check** has been done for this species or genus
  (Onagraceae) - unlike Viscum, which has some precedent context from
  parasitic-plant telomerase-loss literature, this hasn't been searched
  for prior reports at all.

## Reproducing

```bash
# Confirm the TR side is solid (already established elsewhere in this repo):
grep -c "." outputs/tr_tbls/Epilobium_hirsutum.tbl   # non-empty TR hit
grep "^Epilobium_hirsutum" outputs/tr_repeat_correlation/core_template_sets.tsv
grep "^Epilobium_hirsutum" outputs/tr_repeat_correlation/identical_set_array_structure.tsv

# Confirm TERT absence across all 4 standard searches:
for f in outputs/tert_tbls/Epilobium_hirsutum/*.tbl; do grep -vc "^#" "$f"; done   # all 0
wc -l outputs/tert_tblastn/Epilobium_hirsutum.tbl                                  # 0 bytes

# Confirm the congener has a real hit (the positive control):
head -5 outputs/tert_tblastn/Epilobium_ciliatum.tbl

# Run the synteny check:
python3 src/check_gene_synteny_gap.py Epilobium_ciliatum SUPER_16 2595000 2615000 \
  Epilobium_hirsutum outputs/synteny_checks
```
