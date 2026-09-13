# Within-species TR Template variation

**Status: preliminary, dataset still growing.** This is based on a snapshot
of 331 species (out of an eventual ~530+ expected to have multiple loci,
once the current bsub run over the newly-onboarded DToL species finishes -
see `src/pipeline_status.bash`). The headline proportions are unlikely to
shift much with more data, but re-run the query in "Reproducing" below for
the current numbers before citing this anywhere.

## Summary

Yes, there's within-species variation in the TR gene's Template domain, and
it's common, not rare: **95.5% of species with multiple independent
TR-homologous loci (316/331) have at least one locus whose Template
sequence differs from the others.** Most of that variation is exactly what
you'd expect from ordinary duplicated-gene drift (paralogs staying
recognizably canonical while accumulating small differences), but a
substantial minority - **126/331 species (38%)** - show something more
interesting: some loci in the same genome still encode the canonical
repeat while others have already diverged from it. That's consistent with
template pseudogenization caught mid-process, within single genomes.

## Background

This follows directly from the *Carlina vulgaris* investigation
(`notes/carlina_vulgaris_telomere_repeat.md`), which only surfaced because
Carlina happened to have several independent TR-homologous loci and they
disagreed with each other. The rest of the pipeline (`extract_tr_domains.py`
et al.) only ever looks at each species' single best-scoring hit, so it
structurally could not have noticed this pattern in any other species. This
note is the result of checking systematically: `src/extract_tr_template_all_loci.py`
extracts the Template domain from *every* decent TR-homologous locus per
species (E<=1e-4, HMM span>=100 columns - same threshold used dataset-wide),
not just the best one, via real `hmmalign` (not linear-coordinate
interpolation, which would be unreliable for comparing multiple loci against
each other - see the script's docstring).

## Method

For every species with >=2 qualifying loci:

1. Extract each locus (padded genomic region around the hit) and `hmmalign`
   them together against the production HMM.
2. Slice the Template columns (128-174) from the real alignment for each
   locus.
3. Compare the resulting per-locus sequences: are they byte-identical, and
   does each one individually reverse-complement to canonical `TTTAGGG`
   (via the same `find_repeat_in_window` check used elsewhere in the
   pipeline)?

## Results

Of 331 species with >=2 qualifying loci:

| category | species | % |
|---|---|---|
| All loci Template sequences byte-identical | 15 | 4.5% |
| At least one locus differs from the others | 316 | 95.5% |

Breaking the 331 down by canonical-repeat status per locus, rather than raw
sequence identity:

| pattern | species | % |
|---|---|---|
| **All loci canonical** - variation present but every locus still encodes `TTTAGGG` | 200 | 60.4% |
| **Mixed** - some loci canonical, some not, in the same genome | 126 | 38.1% |
| **None canonical** - every locus has diverged | 5 | 1.5% |

The "none canonical" species (*Carlina_vulgaris* is not among them - see
below) are: *Mimulus_peregrinus*, *Fumaria_muralis*,
*Chrysosplenium_oppositifolium*, *Jasione_montana*, *Ranunculus_sardous*.
Three of these (Jasione, Ranunculus, and by implication similar cases) were
already checked in the Carlina work and turned out to be single-nucleotide-
level divergence from canonical rather than a different repeat family - i.e.
"none canonical" by this exact-match test doesn't necessarily mean "no
canonical signal at all," just that it didn't survive a strict test. This
list hasn't been individually re-verified the way Carlina was; treat it as a
lead, not a conclusion.

Notable "mixed" examples (species, loci, canonical/non-canonical split):

| species | loci | canonical | non-canonical |
|---|---|---|---|
| Levisticum_officinale | 47 | 46 | 1 |
| Solanum_nigrum | 40 | 6 | 34 |
| Ajuga_chamaepitys | 34 | 13 | 21 |
| Anagallis_arvensis | 31 | 10 | 21 |
| Platanus_x_hispanica | 22 | 14 | 8 |
| Mentha_aquatica | 15 | 8 | 7 |

*Carlina_vulgaris* itself is "mixed," not "none canonical" as the original
write-up implied - see the corrected `notes/carlina_vulgaris_telomere_repeat.md`
for the update: 6 of its 7 loci carry the AAACTG pattern, but its single
weakest hit (E=9.7e-6, several orders of magnitude less significant than the
other six) reverse-complements cleanly to canonical.

## Interpretation

- The **60.4% "all canonical"** group is the baseline expectation: multi-copy
  ncRNA gene families routinely have several genomic copies (recent
  duplications, retained paralogs) that drift from each other at the
  sequence level without changing which repeat they'd encode.
- The **38.1% "mixed"** group is the interesting one: it's a live snapshot of
  template divergence *within* a genome, at loci that are (by construction)
  close enough in sequence overall to still be found by the same HMM search,
  but have already lost the specific short Template motif's fidelity to the
  ancestral repeat. This is consistent with relaxed selection on
  non-functional/duplicate copies - only one locus needs to stay functional
  (and canonical) for the organism's telomeres to work, so extra copies are
  free to drift, and the Template region (short, and only functionally
  constrained if that particular copy is expressed and incorporated into an
  active telomerase RNP) is a plausible place to see that drift show up
  first.
- The **1.5% "none canonical"** tail is where a *Carlina*-like finding could
  be hiding - but as the Carlina correction above shows, don't trust this
  category at face value without the same depth of follow-up (multiple
  independent loci + genome-wide TIDK positional search) done for Carlina.
  Most of it is probably ordinary divergence past the point where a strict
  matcher can still recognise it, not a repeat-family change.

## Caveats

- **Loci are not independent data points.** Multiple hits in one genome are
  very likely paralogous/duplicated copies of the same ancestral locus
  (or, in some cases, could be allelic/haplotype duplicates in a
  not-fully-collapsed assembly), not independent samples of evolution. A
  species with many loci showing the same divergent pattern is one
  correlated observation, not many.
- **The E<=1e-4 threshold is somewhat arbitrary** and was chosen to match
  the rest of the pipeline for consistency, not tuned specifically for this
  analysis. Loosening or tightening it would change which weak/marginal
  loci count.
- **No haplotype-collapsing has been done.** Some "multiple loci" could be
  the same physical gene appearing on both haplotypes of a phased diploid
  assembly rather than genuinely separate genomic copies - this hasn't been
  checked and could inflate the loci counts for some species (though it
  wouldn't explain the *canonical vs non-canonical* split within a species,
  since the two haplotype copies of the same locus should be nearly
  identical).
- **This snapshot will grow** as the bsub-submitted extraction jobs for the
  newly-onboarded DToL species finish (`src/run_extract_tr_template_all_loci_bsub.bash`) -
  re-run the query below for current numbers.

## Reproducing

```bash
# Already run for all qualifying species (skips species with existing output):
bash src/run_extract_tr_template_all_loci.bash          # local
bash src/run_extract_tr_template_all_loci_bsub.bash      # or via LSF

# The summary numbers above:
python3 - <<'EOF'
import glob, sys
sys.path.insert(0, "src")
from extract_tr_domains import find_repeat_in_window
from collections import defaultdict

species_loci = defaultdict(list)
for f in glob.glob("outputs/tr_template_loci/*.tsv"):
    with open(f) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[1] == "OK" and len(cols) >= 7:
                species_loci[cols[0]].append(cols[6])

all_canonical = mixed = none_canonical = identical = variable = 0
for sp, loci in species_loci.items():
    if len(loci) < 2:
        continue
    if len(set(loci)) == 1:
        identical += 1
    else:
        variable += 1
    hits = [bool(find_repeat_in_window(t, "TTTAGGG")) for t in loci]
    if all(hits): all_canonical += 1
    elif not any(hits): none_canonical += 1
    else: mixed += 1

total = all_canonical + mixed + none_canonical
print(f"multi-locus species: {total}")
print(f"identical: {identical}  variable: {variable}")
print(f"all canonical: {all_canonical}  mixed: {mixed}  none canonical: {none_canonical}")
EOF
```
