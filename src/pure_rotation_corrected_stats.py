#!/usr/bin/env python3
"""Recomputes the headline conservation statistics from
tr_core_template_variation_and_telomeres.md excluding the 4 sets that
turned out to be pure phase-rotations of the modal sequence (0 true
substitutions - tr_core_template_variant_positions.py's is_pure_rotation
flag), which shouldn't count as genuine template divergence. Flagged in
paper_planning.md as needing a check before publication; the correction
turns out to be small (a few percentage points) but real.

Usage:
  pure_rotation_corrected_stats.py <core_template_sets.tsv> <core_template_variant_positions.tsv>
"""
import csv
import sys
from collections import defaultdict


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    sets_rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    vp_rows = list(csv.DictReader(open(sys.argv[2]), delimiter="\t"))

    by_set = defaultdict(list)
    for r in vp_rows:
        by_set[(r["species"], r["core_length"])].append(r)

    pure_rotation_sets = set()
    for key, rs in by_set.items():
        nonmodal = [r for r in rs if r["is_modal"] == "0"]
        if nonmodal and all(r["is_pure_rotation"] == "1" for r in nonmodal):
            pure_rotation_sets.add(key)

    print(f"[pure_rotation_corrected_stats] pure-rotation-only sets "
          f"(reclassified identical): {sorted(pure_rotation_sets)}\n")

    n_total = len(sets_rows)
    n_identical = sum(1 for r in sets_rows if int(r["n_distinct"]) == 1)
    n_identical_corrected = n_identical + len(pure_rotation_sets)
    print(f"All comparable sets: {n_total}")
    print(f"  original:  {n_identical} identical ({100*n_identical/n_total:.1f}%)")
    print(f"  corrected: {n_identical_corrected} identical "
          f"({100*n_identical_corrected/n_total:.1f}%)\n")

    big = [r for r in sets_rows if int(r["n_loci"]) >= 5]
    big_identical = sum(1 for r in big if int(r["n_distinct"]) == 1)
    big_pure_rot = sum(1 for r in big if (r["species"], r["core_length"]) in pure_rotation_sets)
    print(f">=5-loci subset: {len(big)} sets")
    print(f"  original:  {big_identical} identical ({100*big_identical/len(big):.1f}%)")
    print(f"  corrected: {big_identical + big_pure_rot} identical "
          f"({100*(big_identical+big_pure_rot)/len(big):.1f}%)\n")

    differing_species = {r["species"] for r in sets_rows if int(r["n_distinct"]) > 1}
    differing_species_corrected = {
        r["species"] for r in sets_rows
        if int(r["n_distinct"]) > 1 and (r["species"], r["core_length"]) not in pure_rotation_sets
    }
    n_differing_corrected = sum(
        1 for r in sets_rows
        if int(r["n_distinct"]) > 1 and (r["species"], r["core_length"]) not in pure_rotation_sets
    )
    print(f"Differing sets/species:")
    print(f"  original:  44 sets / {len(differing_species)} species")
    print(f"  corrected: {n_differing_corrected} sets / {len(differing_species_corrected)} species")
    print(f"  species dropped entirely: "
          f"{sorted(differing_species - differing_species_corrected)}")


if __name__ == "__main__":
    main()
