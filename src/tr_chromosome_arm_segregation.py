#!/usr/bin/env python3
"""Tests a specific question raised directly by a Zavodnik et al. 2023
(New Phytologist, 10.1111/nph.19110) author in response to this project's
findings: in some species with two TR-supported telomere motifs, the two
variants show chromosome-ARM-specific positions - some arms carry only
one motif, others only the other, with NO co-localisation anywhere in the
genome (their example: Capsicum annuum, their Fig. 4). They flagged this
as awaiting mechanistic clarification - an open question in their own
work, observed qualitatively in a handful of validated species, not
tested systematically.

This reuses the exhaustive per-chromosome-end walk already built for
tr_core_template_array_structure.py (telomere-proper restricted, same
GAP_THRESHOLD, same methodology) but asks a different question of it: not
"is there structure within one array" but "does each chromosome END
(= one arm's telomere) carry the telomere-dominant variant only, the
minor variant(s) only, or both together?" - a direct, dataset-wide,
statistical test of the segregation-vs-colocalisation question, across
all 44 differing sets, rather than a qualitative observation in a
handful of species.

Usage:
  tr_chromosome_arm_segregation.py <core_template_variant_positions.tsv> <positional_dir> [window_bp] [outfile]

Output columns: species, core_length, n_loci, telomere_dominant,
n_arms_checked, n_arms_informative (have >=1 variant present, telomere-
proper), n_dominant_only, n_minor_only, n_both (co-localised), n_neither,
minor_freq (dataset-wide, for context - see
tr_core_template_array_structure.md for the band definitions), pattern
(FULLY_SEGREGATED / PARTIALLY_SEGREGATED / FULLY_COLOCALISED / ABSENT /
TOO_FEW_TO_ASSESS)
"""
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from walk_terminal_repeat_array import prepare_genome, extract_terminal_prepared, tile
from tr_core_template_array_structure import (
    load_sets, strip_label, telomere_proper_prefix,
)
from resolve_weak_sensitivity import all_main_chromosomes

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WINDOW = 5000
MIN_MINOR_FOR_PATTERN = 3  # need at least this many total minor-variant
                            # occurrences dataset-wide before calling a
                            # segregation/colocalisation pattern meaningful


def classify_pattern(n_dominant_only, n_minor_only, n_both, total_minor):
    if total_minor == 0:
        return "ABSENT"
    if total_minor < MIN_MINOR_FOR_PATTERN:
        return "TOO_FEW_TO_ASSESS"
    if n_both == 0 and n_minor_only > 0:
        return "FULLY_SEGREGATED"
    if n_both > 0 and n_minor_only == 0 and n_dominant_only == 0:
        return "FULLY_COLOCALISED"
    return "PARTIALLY_SEGREGATED"


def analyze_species(species, species_sets, positional_dir, window):
    fa, tmpdir = prepare_genome(species)
    results = []
    try:
        for (sp, core_length), s in species_sets:
            modal, variants = s["modal"], s["variants"]
            targets = all_main_chromosomes(sp, modal, positional_dir)

            telomere_counts = Counter()
            per_end_bases = {}
            for chrom, end, approx_length in targets:
                seq, offset, length = extract_terminal_prepared(fa, chrom, end, window, approx_length=approx_length)
                track = tile(seq, variants)
                proper_track = telomere_proper_prefix(track, end)
                bases = set()
                for _, _, label in proper_track:
                    if label:
                        base = strip_label(label)
                        telomere_counts[base] += 1
                        bases.add(base)
                if bases:
                    per_end_bases[(chrom, end)] = bases

            if not telomere_counts:
                results.append(dict(
                    species=sp, core_length=core_length, n_loci=s["n_loci"],
                    telomere_dominant="", n_arms_checked=len(targets), n_arms_informative=0,
                    n_dominant_only=0, n_minor_only=0, n_both=0, n_neither=len(targets),
                    minor_freq=0.0, pattern="ABSENT",
                ))
                continue

            dominant, dom_n = telomere_counts.most_common(1)[0]
            total = sum(telomere_counts.values())
            total_minor = total - dom_n

            n_dominant_only = n_minor_only = n_both = 0
            for (chrom, end), bases in per_end_bases.items():
                has_dom = dominant in bases
                has_minor = len(bases - {dominant}) > 0
                if has_dom and has_minor:
                    n_both += 1
                elif has_dom:
                    n_dominant_only += 1
                elif has_minor:
                    n_minor_only += 1

            n_informative = n_dominant_only + n_minor_only + n_both
            n_neither = len(targets) - n_informative
            minor_freq = total_minor / total if total else 0.0
            pattern = classify_pattern(n_dominant_only, n_minor_only, n_both, total_minor)

            results.append(dict(
                species=sp, core_length=core_length, n_loci=s["n_loci"],
                telomere_dominant=dominant, n_arms_checked=len(targets),
                n_arms_informative=n_informative, n_dominant_only=n_dominant_only,
                n_minor_only=n_minor_only, n_both=n_both, n_neither=n_neither,
                minor_freq=round(minor_freq, 4), pattern=pattern,
            ))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return results


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    variant_positions_path = sys.argv[1]
    positional_dir = sys.argv[2]
    window = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_WINDOW
    outfile = sys.argv[4] if len(sys.argv) > 4 else None

    sets = load_sets(variant_positions_path)
    by_species = defaultdict(list)
    for key, s in sets.items():
        by_species[key[0]].append((key, s))

    all_results = []
    for i, (species, species_sets) in enumerate(sorted(by_species.items()), 1):
        print(f"[tr_chromosome_arm_segregation] ({i}/{len(by_species)}) {species}: "
              f"{len(species_sets)} set(s)", file=sys.stderr)
        try:
            all_results.extend(analyze_species(species, species_sets, positional_dir, window))
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)

    cols = ["species", "core_length", "n_loci", "telomere_dominant", "n_arms_checked",
            "n_arms_informative", "n_dominant_only", "n_minor_only", "n_both", "n_neither",
            "minor_freq", "pattern"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in all_results:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    counts = Counter(r["pattern"] for r in all_results)
    print(f"\n[tr_chromosome_arm_segregation] {len(all_results)} sets analyzed:", file=sys.stderr)
    for k, v in counts.most_common():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
