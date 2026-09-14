#!/usr/bin/env python3
"""Categorize each of the 44 differing sets' telomere_array_structure.tsv
row into an interpretable bucket, combining:
  - minor_freq: what fraction of all telomere-matched tile units belong to
    variant(s) other than the telomere-wide dominant one (NOT the TR-copy
    dominant - see tr_core_template_array_structure.py's docstring for why
    the reference has to be telomere-observed, not gene-copy-count)
  - whether the non-dominant material is a true rare minority or roughly
    co-dominant with the majority (changes what "clustered" means:
    scattered single-unit substitutions vs large sub-array domains)
  - the block/runs evidence already computed

Usage:
  tr_core_template_array_structure_summary.py <telomere_array_structure.tsv> [outfile]
"""
import csv
import sys


def categorize(row):
    n_modal = int(row["n_modal"])
    n_nonmodal = int(row["n_nonmodal"])
    dom_n = int(row["telomere_dom_n"])
    total = n_modal + n_nonmodal
    minor_total = total - dom_n
    minor_freq = minor_total / total if total else 0.0
    max_block = int(row["max_minor_block"])
    runs_n1 = int(row["runs_n1"])

    if minor_total == 0:
        band = "ABSENT"
    elif minor_freq < 0.10:
        band = "RARE"
    elif minor_freq < 0.35:
        band = "MINORITY"
    else:
        band = "CO-DOMINANT/BALANCED"

    if runs_n1 <= 1:
        structure = "TOO_FEW_TO_ASSESS"
    elif max_block >= 2:
        structure = f"BLOCKY (longest run {max_block} units)"
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
            modal=r["modal"], telomere_dominant=r["telomere_dominant"],
            copy_dominant_is_telomere_dominant=r["copy_dominant_is_telomere_dominant"],
            minor_freq=cat["minor_freq"], band=cat["band"], structure=cat["structure"],
            runs_z=r["runs_z"],
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
    print("[tr_core_template_array_structure_summary] structure counts:", Counter(r["structure"] for r in out_rows), file=sys.stderr)


if __name__ == "__main__":
    main()
