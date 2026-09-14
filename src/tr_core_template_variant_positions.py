#!/usr/bin/env python3
"""Full variant breakdown + within-core diff positions for every differing
core-template set (n_distinct > 1 in core_template_sets.tsv).

core_template_sets.tsv only records the top 2 (dominant/minority) variants
per set. This rebuilds the same (species, core_length)-bucketed groups
from all_species.tsv (same anchoring logic as tr_core_template_sets.py)
but keeps the FULL distinct-variant distribution, and for every non-modal
variant reports which position(s) within the core it differs from the
modal (TR-copy-majority) variant.

IMPORTANT: most of these cores are periodic telomere-repeat-family motifs
(e.g. built from an AAC/CCT-like ~3bp unit), so comparing two chunks in
their raw, as-extracted reading frame overstates how different they look
- e.g. AACCCT vs CCCTAA is a naive 5/6-position mismatch but CCCTAA is
exactly AACCCT read starting 2 bases later (0 true substitutions). This
script therefore reports the ROTATION-CORRECTED diff: for each non-modal
variant, it searches all cyclic rotations of the variant against the
fixed modal frame and keeps whichever rotation minimises the Hamming
distance (ties broken toward the smallest |shift|, shift=0 preferred),
alongside the raw (shift=0) distance for transparency. A handful of
pairs turn out to be pure phase reads of the identical sequence (best
diff = 0 at a nonzero shift) - flagged via is_pure_rotation - which are
not point-mutant paralogs at all, just the same repeat unit read from a
different register.

Usage:
  tr_core_template_variant_positions.py <all_species.tsv> [outfile]

Output columns: species, core_length, n_loci, n_distinct, variant,
n_copies, is_modal, target_example, raw_diff_n, best_shift, best_diff_n,
is_pure_rotation, diff_positions (rotation-corrected, comma-separated,
positions in the modal's own frame), diff_bases (e.g. "3:A>T").
"""
import sys
from collections import defaultdict, Counter


def load_rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            yield dict(zip(header, cols))


def build_groups(path):
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
            target=row["target"], chunk=chunk, offset=offset, length=len(chunk),
        ))

    groups = []
    for sp, loci in species_loci.items():
        by_length = defaultdict(list)
        for l in loci:
            by_length[l["length"]].append(l)
        for length, entries in by_length.items():
            if len(entries) < 2:
                continue
            cc = Counter(e["chunk"] for e in entries)
            if len(cc) < 2:
                continue
            groups.append((sp, length, entries, cc))
    return groups


def raw_diff(modal, other):
    return sum(1 for a, b in zip(modal, other) if a != b)


def best_rotation(modal, variant):
    """Search all cyclic rotations of `variant` against the fixed `modal`
    frame; return (best_diff_n, best_shift) minimising Hamming distance,
    ties broken toward the smallest circular |shift| (0 preferred)."""
    n = len(modal)
    order = sorted(range(n), key=lambda k: min(k, n - k))
    best = (n + 1, 0)
    for k in order:
        rot = variant[k:] + variant[:k]
        d = raw_diff(modal, rot)
        if d < best[0]:
            best = (d, k)
    return best


def diff_positions(modal, variant, shift):
    n = len(modal)
    rot = variant[shift:] + variant[:shift]
    positions = [i for i, (a, b) in enumerate(zip(modal, rot)) if a != b]
    diffs = ";".join(f"{i}:{modal[i]}>{rot[i]}" for i in positions)
    return ",".join(str(i) for i in positions), diffs


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else None

    groups = build_groups(path)
    rows = []
    for sp, length, entries, cc in groups:
        modal, _ = cc.most_common(1)[0]
        n_loci = len(entries)
        n_distinct = len(cc)
        for variant, n_copies in cc.most_common():
            target_example = next(e["target"] for e in entries if e["chunk"] == variant)
            is_modal = variant == modal
            if is_modal:
                raw_n = best_diff_n = shift = 0
                is_pure_rot = 0
                dp, db = "", ""
            else:
                raw_n = raw_diff(modal, variant)
                best_diff_n, shift = best_rotation(modal, variant)
                is_pure_rot = int(best_diff_n == 0 and shift != 0)
                dp, db = diff_positions(modal, variant, shift)
            rows.append(dict(
                species=sp, core_length=length, n_loci=n_loci, n_distinct=n_distinct,
                variant=variant, n_copies=n_copies, is_modal=int(is_modal),
                target_example=target_example, raw_diff_n=raw_n, best_shift=shift,
                best_diff_n=best_diff_n, is_pure_rotation=is_pure_rot,
                diff_positions=dp, diff_bases=db,
            ))

    rows.sort(key=lambda r: (-r["n_loci"], r["species"], r["core_length"], -r["n_copies"]))
    cols = ["species", "core_length", "n_loci", "n_distinct", "variant", "n_copies",
            "is_modal", "target_example", "raw_diff_n", "best_shift", "best_diff_n",
            "is_pure_rotation", "diff_positions", "diff_bases"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in rows:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_sets = len({(r["species"], r["core_length"]) for r in rows})
    n_pure_rot = sum(1 for r in rows if r["is_pure_rotation"])
    n_nonmodal = sum(1 for r in rows if not r["is_modal"])
    print(f"[tr_core_template_variant_positions] {n_sets} differing sets, "
          f"{len(rows)} variant rows ({n_nonmodal} non-modal), "
          f"{n_pure_rot} are pure phase-rotations of the modal (0 true substitutions)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
