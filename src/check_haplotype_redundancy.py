#!/usr/bin/env python3
"""Screens all_species.tsv for haplotype-phasing redundancy: some DToL
assemblies are multi-haplotype-phased, and a TR locus can end up
represented more than once under different haplotype-labelled scaffold
names for the SAME true chromosome (e.g. SUPER_12_HAP2 and SUPER_12_HAP4
both carrying a TR hit) - allelic copies of one real locus, not
independent paralogous gene copies. Counting these as independent copies
inflates copy-number statistics.

Handles both haplotype-naming conventions seen across DToL assemblies:
  SUPER_N_HAPx / SUPER_NA_HAPx  (suffix)
  HAPx_SCAFFOLD_N / HAPx_SUPER_N  (prefix)

Usage:
  check_haplotype_redundancy.py <all_species.tsv> [outfile]

Output columns: species, true_chromosome, hap_numbers, n_loci_on_this_chromosome
(one row per true-chromosome that carries a TR locus on more than one
haplotype number - i.e. the flagged/affected rows only).
"""
import csv
import re
import sys
from collections import defaultdict

SUFFIX_RE = re.compile(r'^(.*)_HAP(\d+)$')
PREFIX_RE = re.compile(r'^HAP(\d+)_(.*)$')


def normalize(target):
    m = SUFFIX_RE.match(target)
    if m:
        return m.group(1), m.group(2)
    m = PREFIX_RE.match(target)
    if m:
        return m.group(2), m.group(1)
    return target, None


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else None

    rows = list(csv.DictReader(open(path), delimiter="\t"))
    by_species = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_species[r["species"]][r["target"]].append(r)

    n_any_hap = 0
    flagged = []
    for sp, by_target in by_species.items():
        base_to_haps = defaultdict(set)
        has_hap = False
        for t in by_target:
            base, hap = normalize(t)
            if hap is not None:
                has_hap = True
                base_to_haps[base].add(hap)
        if has_hap:
            n_any_hap += 1
        for base, haps in base_to_haps.items():
            if len(haps) > 1:
                n_loci = sum(len(by_target[t]) for t in by_target
                             if normalize(t)[0] == base)
                flagged.append(dict(species=sp, true_chromosome=base,
                                     hap_numbers=",".join(sorted(haps)),
                                     n_loci_on_this_chromosome=n_loci))

    cols = ["species", "true_chromosome", "hap_numbers", "n_loci_on_this_chromosome"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in flagged:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_species_flagged = len({r["species"] for r in flagged})
    n_total = len(by_species)
    print(f"[check_haplotype_redundancy] {n_any_hap}/{n_total} species have >=1 TR locus "
          f"on a haplotype-labelled scaffold", file=sys.stderr)
    print(f"[check_haplotype_redundancy] {n_species_flagged}/{n_total} species have the "
          f"genuinely inflating pattern (same true chromosome, TR loci on >1 haplotype "
          f"number): {sorted({r['species'] for r in flagged})}", file=sys.stderr)


if __name__ == "__main__":
    main()
