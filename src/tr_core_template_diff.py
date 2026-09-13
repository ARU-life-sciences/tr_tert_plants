#!/usr/bin/env python3
"""Define the actual biological template per TR locus as its position-
anchored best match to a real genome repeat (not the full padded
extraction window from outputs/tr_template_loci, which is deliberately
wide to accommodate length variation across species and includes
flanking scaffold sequence that isn't part of the true minimal
template).

For each species, this:
  1. Takes every locus's matched chunk (from tr_repeat_correlation.py's
     all_species.tsv) and its offset within that locus's extracted
     Template window (template.find(chunk)) - the chunk's genuine
     anchored position, not just its content.
  2. Compares pairs of same-species loci ONLY when their chunks are the
     same length AND anchor at the same (or near-same) offset - i.e.
     genuinely the same structural position in the aligned Template
     domain, not a coincidental same-length match from unrelated
     positions (the mistake that produced a misleading comparison
     earlier in this analysis).
  3. Reports the true core-template edit distance for those pairs, so
     "do two TR gene copies differ at the templating position" can be
     answered directly rather than inferred from whole-window Hamming
     distance (which conflates templating-region and flanking-scaffold
     divergence).

Usage:
  tr_core_template_diff.py <all_species.tsv> [max_offset_diff] [outfile]
"""
import sys
from collections import defaultdict
from itertools import combinations

DEFAULT_MAX_OFFSET_DIFF = 1  # allow the anchor to differ by at most this many columns


def load_rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            yield dict(zip(header, cols))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    max_offset_diff = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_MAX_OFFSET_DIFF
    outfile = sys.argv[3] if len(sys.argv) > 3 else None

    species_loci = defaultdict(list)
    for row in load_rows(path):
        mlen = int(row["match_len"])
        if mlen == 0:
            continue
        template = row["template"]
        chunk = row["match_chunk"]
        offset = template.find(chunk) if row["match_orientation"] == "direct" else None
        if offset is None:
            continue  # skip revcomp-orientation rows - offset isn't meaningful against `template` directly
        species_loci[row["species"]].append(dict(
            target=row["target"], coords=row["coords"], evalue=row["evalue"],
            template=template, chunk=chunk, offset=offset, p=float(row["shuffle_p"]),
        ))

    out_rows = []
    for sp, loci in species_loci.items():
        for l1, l2 in combinations(loci, 2):
            if l1["target"] == l2["target"] and l1["coords"] == l2["coords"]:
                continue
            if len(l1["chunk"]) != len(l2["chunk"]):
                continue
            if abs(l1["offset"] - l2["offset"]) > max_offset_diff:
                continue  # not anchored at the same structural position - not comparable
            d = sum(x != y for x, y in zip(l1["chunk"], l2["chunk"]))
            out_rows.append((sp, l1["target"], l1["coords"], l1["evalue"], l1["chunk"], l1["offset"],
                              l2["target"], l2["coords"], l2["evalue"], l2["chunk"], l2["offset"], d))

    out_rows.sort(key=lambda r: r[11])  # smallest core edit distance first

    cols = ["species", "locus1_target", "locus1_coords", "locus1_evalue", "locus1_core",
            "locus1_offset", "locus2_target", "locus2_coords", "locus2_evalue", "locus2_core",
            "locus2_offset", "core_edit_distance"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in out_rows:
        print("\t".join(str(x) for x in r), file=out)
    if outfile:
        out.close()
        print(f"[tr_core_template_diff] {len(out_rows)} position-anchored comparable pairs -> {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
