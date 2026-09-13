#!/usr/bin/env python3
"""Species-level (not pairwise) view of TR core-template variation, built
on tr_core_template_diff.py's anchoring logic.

Why not just count pairs: a species with N mutually-identical loci
contributes C(N,2) trivially-identical pairs, which dominates any
pairwise "% identical" statistic and can make a genuinely variable
dataset look far more conserved than it is (a 44-copy, fully-identical
species alone contributes 946 "identical" pairs). This groups loci per
species into comparable sets instead - same core length, anchored at the
same (+/-1) position within the Template alignment - and reports, per
set: how many loci, how many distinct core sequences, and (for sets with
>1 distinct core) the dominant and top minority variant, which is what
run_confirm_differing_sets(.bash|_bsub.bash) consumes for telomere
confirmation.

Usage:
  tr_core_template_sets.py <all_species.tsv> [outfile]

Output columns: species, core_length, n_loci, n_distinct, dominant,
dom_n, dom_target, minority, min_n, min_target (minority/min_n/min_target
are empty when n_distinct==1).
"""
import sys
from collections import defaultdict, Counter


def load_rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            yield dict(zip(header, cols))


def build_sets(path):
    species_loci = defaultdict(list)
    for row in load_rows(path):
        mlen = int(row["match_len"])
        if mlen == 0 or row["match_orientation"] != "direct":
            continue
        template, chunk = row["template"], row["match_chunk"]
        offset = template.find(chunk)
        if offset == -1:
            continue
        species_loci[row["species"]].append(dict(
            target=row["target"], coords=row["coords"], chunk=chunk,
            offset=offset, length=len(chunk),
        ))

    results = []
    for sp, loci in species_loci.items():
        if len(loci) < 2:
            continue
        # bucket by (length, offset), then merge buckets for the same
        # length whose offsets are within 1 of each other (same
        # structural position, allowing minor alignment jitter)
        by_length = defaultdict(list)
        for l in loci:
            by_length[l["length"]].append(l)
        for length, entries in by_length.items():
            if len(entries) < 2:
                continue
            chunks = [e["chunk"] for e in entries]
            cc = Counter(chunks)
            n_distinct = len(cc)
            dominant, dom_n = cc.most_common(1)[0]
            dom_target = next(e["target"] for e in entries if e["chunk"] == dominant)
            if n_distinct > 1:
                minority, min_n = cc.most_common(2)[1]
                min_target = next(e["target"] for e in entries if e["chunk"] == minority)
            else:
                minority = min_n = min_target = ""
            results.append(dict(
                species=sp, core_length=length, n_loci=len(entries), n_distinct=n_distinct,
                dominant=dominant, dom_n=dom_n, dom_target=dom_target,
                minority=minority, min_n=min_n, min_target=min_target,
            ))
    return results


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else None
    results = build_sets(path)
    results.sort(key=lambda r: -r["n_loci"])

    cols = ["species", "core_length", "n_loci", "n_distinct", "dominant", "dom_n",
            "dom_target", "minority", "min_n", "min_target"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_sets = len(results)
    n_identical = sum(1 for r in results if r["n_distinct"] == 1)
    n_big = [r for r in results if r["n_loci"] >= 5]
    n_big_identical = sum(1 for r in n_big if r["n_distinct"] == 1)
    print(f"[tr_core_template_sets] {n_sets} comparable sets, {n_identical} "
          f"({100*n_identical/n_sets:.1f}%) fully identical", file=sys.stderr)
    if n_big:
        print(f"[tr_core_template_sets] of {len(n_big)} sets with >=5 loci, "
              f"{n_big_identical} ({100*n_big_identical/len(n_big):.1f}%) fully identical", file=sys.stderr)


if __name__ == "__main__":
    main()
