#!/usr/bin/env python3
"""Exhaustively resolve every WEAK_SENSITIVITY set from
summarize_repeat_confirmations.py - cases where the dominant TR-copy
variant is confirmed telomeric by confirm_repeat_positional.py but the
minority variant shows no signal by that coarse (window-count-threshold)
method. As shown for Cornus_sanguinea and Arctium_lappa, a real,
low-frequency, interspersed minority variant can fail that threshold
while still being genuinely present - only base-pair-resolved sequence
walking (walk_terminal_repeat_array.py) can settle it.

Checks EVERY main chromosome (both 5' and 3' ends - not just the
top-N strongest for the dominant variant): a minority variant can be
present at only a subset of chromosomes unrelated to which ones are
strongest for the dominant one, so any partial-coverage check risks a
false GENUINELY_ABSENT call (confirmed: an earlier top-3-only pass
missed real presence in 2/6 of the cases that turned out positive).
One genome decompression per species (prepare_genome), reused across
every chromosome-end checked for that species - not one per fetch.

Usage:
  resolve_weak_sensitivity.py <differing_set_confirmations.tsv> [positional_dir] [window_bp] [outfile]

Requires the dominant variant's raw scan
(<positional_dir>/<species>_<dominant>_telomeric_repeat_windows.tsv,
from confirm_repeat_positional.py) to already exist, to enumerate that
species' main chromosomes and give an approximate length per chromosome
(for a fast 3'-end fetch - see walk_terminal_repeat_array.extract_terminal_prepared).
"""
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from walk_terminal_repeat_array import prepare_genome, extract_terminal_prepared, tile

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WINDOW = 5000
MIN_MAIN_WINDOWS = 100
MAX_MAIN_CHROMS = 40
TIDK_WINDOW_BP = 5000  # matches confirm_repeat_positional.py's default -w


def all_main_chromosomes(species, dominant_variant, positional_dir):
    """Every main chromosome (both ends) for this species, with an
    approximate length per chromosome (n_windows * TIDK_WINDOW_BP) for
    the fast 3'-end fetch."""
    path = Path(positional_dir) / f"{species}_{dominant_variant}_telomeric_repeat_windows.tsv"
    if not path.exists():
        return []
    data = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            data[row["id"]].append(int(row["window"]))
    sized = sorted(data.items(), key=lambda kv: -len(kv[1]))
    main = [(c, w) for c, w in sized if len(w) >= MIN_MAIN_WINDOWS][:MAX_MAIN_CHROMS]

    out = []
    for chrom, windows in main:
        approx_length = max(windows) if windows else len(windows) * TIDK_WINDOW_BP
        out.append((chrom, "5prime", None))
        out.append((chrom, "3prime", approx_length))
    return out


def resolve(species, dominant, minority, positional_dir, window):
    targets = all_main_chromosomes(species, dominant, positional_dir)
    if not targets:
        return dict(species=species, dominant=dominant, minority=minority,
                    verdict="NO_DOMINANT_SCAN", n_chroms_checked=0, dom_n=0, min_n=0,
                    min_freq=0.0, found_at="")

    fa, tmpdir = prepare_genome(species)
    try:
        total_dom, total_min, found_at = 0, 0, []
        n_checked_chroms = len({c for c, _, _ in targets})
        for chrom, end, approx_length in targets:
            seq, offset, length = extract_terminal_prepared(fa, chrom, end, window, approx_length=approx_length)
            track = tile(seq, [dominant, minority])
            counts = {}
            for _, _, label in track:
                if label:
                    counts[label] = counts.get(label, 0) + 1
            dom_n = counts.get(dominant, 0) + counts.get(f"revcomp({dominant})", 0)
            min_n = counts.get(minority, 0) + counts.get(f"revcomp({minority})", 0)
            total_dom += dom_n
            total_min += min_n
            if min_n > 0:
                found_at.append(f"{chrom}:{end}({min_n})")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    total = total_dom + total_min
    freq = total_min / total if total else 0.0
    verdict = "CONFIRMED_PRESENT" if total_min > 0 else "GENUINELY_ABSENT"
    return dict(species=species, dominant=dominant, minority=minority, verdict=verdict,
                n_chroms_checked=n_checked_chroms, dom_n=total_dom, min_n=total_min,
                min_freq=round(freq, 4), found_at=";".join(found_at))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    confirmations_path = sys.argv[1]
    positional_dir = sys.argv[2] if len(sys.argv) > 2 else str(PROJECT_ROOT / "outputs" / "tr_repeat_correlation" / "positional")
    window = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_WINDOW
    outfile = sys.argv[4] if len(sys.argv) > 4 else None

    rows = []
    with open(confirmations_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["classification"] != "WEAK_SENSITIVITY":
                continue
            print(f"[resolve_weak_sensitivity] {r['species']} {r['dominant']} vs {r['minority']} "
                  f"(exhaustive: all main chromosomes, both ends)", file=sys.stderr)
            result = resolve(r["species"], r["dominant"], r["minority"], positional_dir, window)
            rows.append(result)
            print(f"  -> {result['verdict']}, checked {result['n_chroms_checked']} chromosomes, "
                  f"min_n={result['min_n']}", file=sys.stderr)

    cols = ["species", "dominant", "minority", "verdict", "n_chroms_checked", "dom_n", "min_n", "min_freq", "found_at"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in rows:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_confirmed = sum(1 for r in rows if r["verdict"] == "CONFIRMED_PRESENT")
    print(f"\n[resolve_weak_sensitivity] {n_confirmed}/{len(rows)} WEAK_SENSITIVITY sets "
          f"confirmed present by exhaustive direct sequence walking", file=sys.stderr)


if __name__ == "__main__":
    main()
