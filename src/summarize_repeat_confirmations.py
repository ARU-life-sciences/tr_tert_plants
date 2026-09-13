#!/usr/bin/env python3
"""Summarize confirm_repeat_positional.py output for every differing set
from tr_core_template_sets.py: does the minority TR-gene-copy variant
show real chromosome-terminal enrichment, and does it overlap the same
chromosomes as the dominant variant?

Re-implements confirm_repeat_positional.py's own enrichment call (same
thresholds) directly against its raw per-chromosome TSVs, so this stays
consistent with that script without re-running tidk.

Classifies each set:
  STRONG_OVERLAP   both variants enriched, and they substantially share
                    the same chromosome set (overlap >= 80% of the
                    smaller side)
  PARTIAL_OVERLAP   both enriched, but less shared overlap
  WEAK_SENSITIVITY  dominant enriched, minority not (or only marginally) -
                    per notes/tr_core_template_variation_and_telomeres.md,
                    this can still be a real positive that the window-count
                    threshold is too coarse to see (see e.g. Cornus_sanguinea,
                    confirmed present by direct base-pair walking despite
                    scoring here as weak) - treat as "needs the walk tool",
                    not as a confirmed negative
  LABEL_FLIP        dominant NOT enriched but minority IS - the by-TR-copy-
                    count "dominant" label doesn't match which variant is
                    actually the genome's real telomere repeat
  NEITHER           neither variant shows real signal - most likely a
                    genuine negative (see Ajuga_chamaepitys's AACCTAATC)
  MISSING_DATA      one or both raw scans not found - run
                    run_confirm_differing_sets(.bash|_bsub.bash) first

Usage:
  summarize_repeat_confirmations.py <core_template_sets.tsv> [positional_dir] [outfile]
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MIN_TERMINAL_COUNT = 10
MIN_FOLD_OVER_MEDIAN = 5
MIN_MAIN_WINDOWS = 100
MAX_MAIN_CHROMS = 40


def enriched_chromosomes(species, variant, positional_dir):
    path = Path(positional_dir) / f"{species}_{variant}_telomeric_repeat_windows.tsv"
    if not path.exists():
        return None
    data = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            data[row["id"]].append((int(row["window"]), int(row["forward_repeat_number"]), int(row["reverse_repeat_number"])))
    sized = sorted(data.items(), key=lambda kv: -len(kv[1]))
    main_ids = {c for c, w in sized if len(w) >= MIN_MAIN_WINDOWS}
    if len(main_ids) > MAX_MAIN_CHROMS:
        main_ids = {c for c, _ in sized[:MAX_MAIN_CHROMS]}

    enriched = set()
    for chrom, windows in data.items():
        if chrom not in main_ids:
            continue
        windows = sorted(windows)
        totals = [f + r for _, f, r in windows]
        nonzero = [t for t in totals if t > 0]
        median = sorted(nonzero)[len(nonzero) // 2] if nonzero else 0

        def enr(c):
            return c >= MIN_TERMINAL_COUNT and (median == 0 or c >= MIN_FOLD_OVER_MEDIAN * median)

        if enr(totals[0]) or enr(totals[-1]):
            enriched.add(chrom)
    return enriched


def classify(dom_chroms, min_chroms):
    if not dom_chroms and not min_chroms:
        return "NEITHER"
    if not dom_chroms and min_chroms:
        return "LABEL_FLIP"
    if dom_chroms and not min_chroms:
        return "WEAK_SENSITIVITY"
    overlap = dom_chroms & min_chroms
    frac = len(overlap) / min(len(dom_chroms), len(min_chroms))
    return "STRONG_OVERLAP" if frac >= 0.8 else "PARTIAL_OVERLAP"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sets_path = sys.argv[1]
    positional_dir = sys.argv[2] if len(sys.argv) > 2 else str(PROJECT_ROOT / "outputs" / "tr_repeat_correlation" / "positional")
    outfile = sys.argv[3] if len(sys.argv) > 3 else None

    rows = []
    with open(sets_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if int(r["n_distinct"]) <= 1:
                continue
            dom = enriched_chromosomes(r["species"], r["dominant"], positional_dir)
            mino = enriched_chromosomes(r["species"], r["minority"], positional_dir)
            if dom is None or mino is None:
                rows.append((r["species"], r["dominant"], r["minority"], "MISSING_DATA", "", "", ""))
                continue
            rows.append((r["species"], r["dominant"], r["minority"],
                         classify(dom, mino), len(dom), len(mino), len(dom & mino)))

    cols = ["species", "dominant", "minority", "classification", "dom_chroms", "min_chroms", "overlap_chroms"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in rows:
        print("\t".join(str(x) for x in r), file=out)
    if outfile:
        out.close()

    from collections import Counter
    counts = Counter(r[3] for r in rows)
    print(f"[summarize_repeat_confirmations] {len(rows)} differing sets classified:", file=sys.stderr)
    for k, v in counts.most_common():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
