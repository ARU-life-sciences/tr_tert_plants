# Telomere-independent Template prediction

Can the TR gene's core Template be located without reference to the
species' own telomere repeat (i.e. without TIDK)? Checked against the
1073 telomere-verified cores in `outputs/tr_repeat_correlation/all_species.tsv`
(significant, direct-orientation matches only):

- The **G-rich-5'-end -> Template boundary** (HMM column 128, found by
  pure cross-species conservation - no repeat data used) is a tight
  anchor: the real core starts within 3bp of it in 91.1% of loci.
- The **Template -> Conserved-region boundary** is *not* useful as an
  end-anchor: the real core ends a median of 31bp before it.
- Core length is fairly consistent: median 12bp (IQR 11-13, p10-p90
  9-14).

**Rule**: start 2bp after the G-rich boundary, take 12bp. Mean
intersection-over-union against the true core: **0.836** across all 1073
loci, using no telomere data at all. Implemented in
`src/predict_tr_template.py`.

Caveat: this is an average, not a per-locus guarantee - true cores range
7-17bp, and accuracy for genuinely divergent species (e.g. Carlina) is
untested. Good enough for a first-pass call where TIDK data isn't
available; not a replacement for the telomere-matched method where it is.
