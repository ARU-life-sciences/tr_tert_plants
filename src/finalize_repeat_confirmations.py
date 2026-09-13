#!/usr/bin/env python3
"""Merge summarize_repeat_confirmations.py's coarse-method classification
with resolve_weak_sensitivity.py's exhaustive walk-based resolution of
the WEAK_SENSITIVITY cases, into one final table.

Keeps the coarse classification column as-is (it documents what the
window-count-threshold method alone finds, which matters
methodologically - see notes/tr_core_template_variation_and_telomeres.md's
calibration discussion), and adds a `final_verdict` column:
  - STRONG_OVERLAP / PARTIAL_OVERLAP -> CONFIRMED_PRESENT (both variants
    already directly observed by the coarse method)
  - WEAK_SENSITIVITY -> the resolve_weak_sensitivity.py verdict
    (CONFIRMED_PRESENT or GENUINELY_ABSENT, from exhaustively walking
    every main chromosome's both ends)
  - LABEL_FLIP -> LABEL_FLIP (a different, unresolved phenomenon - not
    settled by this exhaustive check, left as-is)

Usage:
  finalize_repeat_confirmations.py <differing_set_confirmations.tsv> <weak_sensitivity_resolved.tsv> [outfile]
"""
import sys


def load_rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    confirmations_path, resolved_path = sys.argv[1], sys.argv[2]
    outfile = sys.argv[3] if len(sys.argv) > 3 else None

    resolved = {(r["species"], r["dominant"], r["minority"]): r for r in load_rows(resolved_path)}
    confirmations = load_rows(confirmations_path)

    out_rows = []
    for r in confirmations:
        key = (r["species"], r["dominant"], r["minority"])
        if r["classification"] in ("STRONG_OVERLAP", "PARTIAL_OVERLAP"):
            final_verdict = "CONFIRMED_PRESENT"
            min_freq = ""
        elif r["classification"] == "WEAK_SENSITIVITY":
            res = resolved.get(key)
            final_verdict = res["verdict"] if res else "UNRESOLVED"
            min_freq = res["min_freq"] if res else ""
        else:
            final_verdict = r["classification"]  # LABEL_FLIP, NEITHER, MISSING_DATA
            min_freq = ""
        out_rows.append(dict(r, final_verdict=final_verdict, min_freq=min_freq))

    cols = ["species", "dominant", "minority", "classification", "final_verdict", "min_freq",
            "dom_chroms", "min_chroms", "overlap_chroms"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in out_rows:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    from collections import Counter
    counts = Counter(r["final_verdict"] for r in out_rows)
    print(f"[finalize_repeat_confirmations] {len(out_rows)} sets, final verdicts:", file=sys.stderr)
    for k, v in counts.most_common():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
