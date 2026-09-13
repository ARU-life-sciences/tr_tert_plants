#!/usr/bin/env python3
"""Generalization of the positional confirmation step used on Carlina
vulgaris (notes/carlina_vulgaris_telomere_repeat.md): run `tidk search`
for a specific repeat string directly on a species' genome and check
whether it shows a genuine chromosome-terminal telomere signature (counts
spiking at the very first/last windows of chromosomes, decaying within a
few windows moving inward) rather than being merely abundant
genome-wide - tidk explore's ranked list only tells you a repeat is
common within its --distance-restricted window, not that it's
specifically terminal, so this is the actual positional check.

Usage:
  confirm_repeat_positional.py <species> <repeat_string> [outdir] [window_bp]

Writes <outdir>/<species>_<repeat>_telomeric_repeat_windows.tsv (raw tidk
output) and prints a per-chromosome summary + a verdict.
"""
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TIDK = "tidk"
DEFAULT_WINDOW = 5000
# A chromosome counts as showing terminal enrichment for this repeat if
# its first or last window's count is at least this many multiples of the
# chromosome's own median non-zero window count, AND at least this many
# raw counts (avoids calling noise on a chromosome with ~0 signal
# everywhere "enriched" just because 1 window has a count of 2).
MIN_TERMINAL_COUNT = 10
MIN_FOLD_OVER_MEDIAN = 5


def load_genome_path(species):
    with open(PROJECT_ROOT / "inputs" / "dtol_plant_paths.txt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == species:
                return parts[1]
    return None


def run_tidk_search(species, repeat, genome_path, outdir, window):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    label = f"{species}_{repeat}"
    subprocess.run(
        [TIDK, "search", "-s", repeat, "-w", str(window), "-o", label,
         "-d", str(outdir), "--log", genome_path],
        check=True, capture_output=True,
    )
    return outdir / f"{label}_telomeric_repeat_windows.tsv"


def analyse(tsv_path):
    data = defaultdict(list)
    with open(tsv_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            data[row["id"]].append((int(row["window"]),
                                     int(row["forward_repeat_number"]),
                                     int(row["reverse_repeat_number"])))

    # "Main" chromosomes are identified by scale, not naming convention -
    # not every assembly uses the "SUPER_N" convention (e.g. some only
    # reach scaffold-level), so take sequences with at least MIN_WINDOWS
    # windows (i.e. >= MIN_WINDOWS * window_bp long) rather than pattern-
    # matching names, up to a generous cap so a handful of huge
    # unplaced scaffolds in a fragmented assembly don't get treated as
    # chromosomes.
    MIN_WINDOWS = 100  # >=500kb at the default 5kb window
    MAX_MAIN_CHROMS = 40
    sized = sorted(data.items(), key=lambda kv: -len(kv[1]))
    main_ids = {chrom for chrom, windows in sized
                if len(windows) >= MIN_WINDOWS}
    if len(main_ids) > MAX_MAIN_CHROMS:
        main_ids = {chrom for chrom, _ in sized[:MAX_MAIN_CHROMS]}

    results = {}
    for chrom, windows in data.items():
        if chrom not in main_ids:
            continue
        windows = sorted(windows)
        totals = [f + r for _, f, r in windows]
        nonzero = [t for t in totals if t > 0]
        median = sorted(nonzero)[len(nonzero) // 2] if nonzero else 0
        first_total = totals[0]
        last_total = totals[-1]
        first_fwd, first_rev = windows[0][1], windows[0][2]
        last_fwd, last_rev = windows[-1][1], windows[-1][2]

        def enriched(count):
            return count >= MIN_TERMINAL_COUNT and (median == 0 or count >= MIN_FOLD_OVER_MEDIAN * median)

        five_prime = enriched(first_total)
        three_prime = enriched(last_total)
        results[chrom] = dict(first_fwd=first_fwd, first_rev=first_rev,
                               last_fwd=last_fwd, last_rev=last_rev,
                               median_window=median,
                               five_prime_enriched=five_prime,
                               three_prime_enriched=three_prime)
    return results


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    species = sys.argv[1]
    repeat = sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else "outputs/tr_repeat_correlation/positional"
    window = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_WINDOW

    genome_path = load_genome_path(species)
    if genome_path is None:
        sys.exit(f"species {species!r} not found in inputs/dtol_plant_paths.txt")

    print(f"[confirm_repeat_positional] {species}: running tidk search -s {repeat}", file=sys.stderr)
    tsv_path = run_tidk_search(species, repeat, genome_path, outdir, window)
    results = analyse(tsv_path)

    n_chroms = len(results)
    n_enriched = sum(1 for r in results.values() if r["five_prime_enriched"] or r["three_prime_enriched"])

    print(f"species\trepeat\tchromosome\tfirst_fwd\tfirst_rev\tlast_fwd\tlast_rev\tmedian_window\t5prime_enriched\t3prime_enriched")
    for chrom, r in sorted(results.items()):
        print(f"{species}\t{repeat}\t{chrom}\t{r['first_fwd']}\t{r['first_rev']}\t"
              f"{r['last_fwd']}\t{r['last_rev']}\t{r['median_window']}\t"
              f"{r['five_prime_enriched']}\t{r['three_prime_enriched']}")

    verdict = "TERMINAL_ENRICHED" if n_enriched >= 2 else ("WEAK" if n_enriched == 1 else "NOT_TERMINAL")
    print(f"[confirm_repeat_positional] {species} / {repeat}: {n_enriched}/{n_chroms} chromosomes show "
          f"terminal enrichment -> {verdict}", file=sys.stderr)


if __name__ == "__main__":
    main()
