#!/usr/bin/env python3
"""Categorize each telomere_array_structure.tsv row (or the identical-set
control group's) into an interpretable bucket, combining:
  - minor_freq: what fraction of all telomere-proper-matched tile units
    belong to variant(s) other than the telomere-dominant one (NOT the
    TR-copy dominant - see tr_core_template_array_structure.py's
    docstring for why the reference has to be telomere-observed)
  - whether the non-dominant material is a true rare minority or roughly
    co-dominant with the majority
  - the structural classification, from the formal runs-test statistic
    FIRST (distinguishing REGULAR_ALTERNATION - strongly positive z, a
    candidate HOR signature - from ordinary CLUSTERED - strongly negative
    z or, when the z-test is unpowered, a directly-observed block of >=2)
    rather than defaulting to "any block found" before checking the
    statistic, which would silently absorb regular-alternation cases into
    a generic "clustered" label.

Uses the telomere-PROPER-restricted columns (proper_*) throughout, not
the full-window ones - see tr_core_template_array_structure.py's
telomere_proper_prefix() for why the full window can include
subtelomeric/interstitial sequence beyond the true terminal array.

Usage:
  tr_core_template_array_structure_summary.py <telomere_array_structure.tsv> [outfile]
"""
import csv
import sys


def categorize(row):
    n_modal = int(row["n_modal"]) if "n_modal" in row else 0
    dom_n = int(row["proper_dom_n"])
    minor_total = int(row["proper_n_minor"])
    total = dom_n + minor_total
    minor_freq = minor_total / total if total else 0.0
    max_block = int(row["proper_max_block"])
    z = row["proper_runs_z"]
    z = float(z) if z not in ("None", "") else None

    if minor_total == 0:
        band = "ABSENT"
    elif minor_freq < 0.10:
        band = "RARE"
    elif minor_freq < 0.35:
        band = "MINORITY"
    else:
        band = "CO-DOMINANT/BALANCED"

    if minor_total == 0:
        structure = "ABSENT"
    elif z is not None and z > 1.96:
        structure = f"REGULAR_ALTERNATION (block cap {max_block}, z={z:+.1f})"
    elif z is not None and z < -1.96:
        structure = f"CLUSTERED (runs test, block {max_block}, z={z:+.1f})"
    elif max_block >= 2:
        structure = f"BLOCKY (longest run {max_block} units)"
    elif int(row["runs_n1"]) <= 1 if "runs_n1" in row else False:
        structure = "TOO_FEW_TO_ASSESS"
    else:
        structure = "ISOLATED (never adjacent to itself)"

    return dict(minor_freq=round(minor_freq, 4), band=band, structure=structure)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else None

    rows = list(csv.DictReader(open(path), delimiter="\t"))
    out_rows = []
    for r in rows:
        cat = categorize(r)
        out_rows.append(dict(
            species=r["species"], core_length=r["core_length"], n_loci=r["n_loci"],
            modal=r["modal"], telomere_dominant=r["proper_dominant"],
            copy_dominant_is_telomere_dominant=r["proper_copy_is_dominant"],
            minor_freq=cat["minor_freq"], band=cat["band"], structure=cat["structure"],
            runs_z=r["proper_runs_z"],
        ))
    out_rows.sort(key=lambda r: r["minor_freq"])

    cols = ["species", "core_length", "n_loci", "modal", "telomere_dominant",
            "copy_dominant_is_telomere_dominant", "minor_freq", "band", "structure", "runs_z"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in out_rows:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    from collections import Counter
    print("[tr_core_template_array_structure_summary] band counts:", Counter(r["band"] for r in out_rows), file=sys.stderr)
    struct_kind = Counter(r["structure"].split(" (")[0] for r in out_rows)
    print("[tr_core_template_array_structure_summary] structure counts:", struct_kind, file=sys.stderr)


if __name__ == "__main__":
    main()
