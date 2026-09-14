#!/usr/bin/env python3
"""Formalizes the per-locus shuffle-control significance test in
all_species.tsv with a proper multiple-testing correction (Benjamini-
Hochberg FDR), rather than the informal "raw p==0.0 (0/200 shuffles)"
threshold used elsewhere. This was flagged in paper_planning.md as
needed before any significance claim from the ~1450-locus screen goes in
a paper.

Empirical shuffle p-values of exactly 0.0 (0/200 shuffles reached the
real match length) are floored to 1/(n_shuffles+1) before correction -
the finest resolution 200 shuffles can actually distinguish; a true
p-value could be smaller, but this is unmeasurable at this shuffle count
(a real caveat, reported below - q < 0.001 is categorically unreachable
with only 200 shuffles per locus).

Usage:
  multiple_testing_correction.py <all_species.tsv>
"""
import csv
import sys

import numpy as np
from scipy.stats import false_discovery_control


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    matched = [r for r in rows if int(r["match_len"]) > 0]

    n_shuffles = int(matched[0]["n_shuffles"])
    floor_p = 1 / (n_shuffles + 1)

    pvals = np.array([floor_p if float(r["shuffle_p"]) == 0.0 else float(r["shuffle_p"])
                       for r in matched])
    qvals = false_discovery_control(pvals, method="bh")

    n_naive = sum(1 for p in pvals if p == floor_p)
    print(f"[multiple_testing_correction] {len(rows)} total loci, {len(matched)} with a "
          f"match (tested)")
    print(f"[multiple_testing_correction] shuffle floor p-value ({n_shuffles} shuffles): "
          f"{floor_p:.5f} - q < {floor_p:.3f} is unreachable at this shuffle count")
    print(f"[multiple_testing_correction] naive uncorrected (raw p=0.0 only): "
          f"{n_naive}/{len(matched)} ({100*n_naive/len(matched):.1f}%)")
    for alpha in (0.05, 0.01, 0.001):
        n_sig = int((qvals < alpha).sum())
        print(f"[multiple_testing_correction] BH-FDR q < {alpha}: {n_sig}/{len(matched)} "
              f"loci significant ({100*n_sig/len(matched):.1f}%)")


if __name__ == "__main__":
    main()
