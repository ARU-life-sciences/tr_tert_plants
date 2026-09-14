#!/usr/bin/env python3
"""For species with FULLY IDENTICAL TR gene copies (n_distinct==1 in
core_template_sets.tsv - no paralog sequence divergence at all), test
whether their genome's own tidk-discovered repeat candidates still show
a plausible second, close variant - completely independent of TR gene
content. This is the control group for asking whether blocky/alternating
telomere array structure (as found in the 44 TR-paralog-differing sets)
can also arise WITHOUT any TR paralog divergence, i.e. a direct test of
Belyayev et al. 2023's DNA-recombination-only mechanism against our
TR-paralog-linked framing.

For each identical set: finds the dominant TR-derived sequence's match
among that species' outputs/tidk/<species>.tidk.tsv candidates (sanity
check - it should be there, since match_chunk was originally selected by
correlating against this same tidk list), then searches the REST of the
tidk candidate list for the best same-length, low-edit-distance (after
cyclic rotation correction, same logic as tr_core_template_variant_positions.py)
"close second candidate" with a non-trivial count - a candidate "derivative
monomer" in Belyayev's sense, discovered independent of any TR paralog.

Usage:
  find_identical_set_second_candidates.py <core_template_sets.tsv> <tidk_dir> [outfile] [variant_positions_outfile]

Output columns: species, core_length, n_loci, dominant, dominant_tidk_count,
second, second_tidk_count, best_diff_n, best_shift
(rows only for sets where a viable second candidate was found; pure
rotations of the dominant - best_diff_n==0 - are excluded upstream, not
flagged, since they wouldn't be a real second variant at all)

If variant_positions_outfile is given, also writes a
tr_core_template_array_structure.py-compatible variant-positions TSV (two
rows per set: modal=dominant, and the second candidate) directly usable as
input to tr_identical_set_array_structure.py, so the control-group pipeline
doesn't need a separate manual adapter step.
"""
import csv
import sys
from pathlib import Path

MAX_DIFF = 3  # tightened from an initial pass at 5, which let through obvious noise
              # (near-homopolymers, unrelated satellites within Hamming-ball-by-chance)
MIN_SECOND_COUNT = 200  # tidk count_repeat_runs_gt_100 floor - not just noise
# NOTE: this is a coarse pre-filter only - tidk explore's candidate list is
# terminal-proximal-restricted (--distance 0.01) but genome-wide across all
# chromosome ends, so it can still include real-but-unrelated
# terminal/subtelomeric repeats (rDNA, satellites) that happen to pass this
# filter by chance. The decisive filter is the actual base-pair-resolved
# terminal walk in the next stage - candidates that don't show real presence
# there are rejected, exactly as for the WEAK_SENSITIVITY/LABEL_FLIP cases.


def best_rotation(a, b):
    n = len(a)
    order = sorted(range(n), key=lambda k: min(k, n - k))
    best = (n + 1, 0)
    for k in order:
        rot = b[k:] + b[:k]
        d = sum(1 for x, y in zip(a, rot) if x != y)
        if d < best[0]:
            best = (d, k)
    return best


def load_tidk(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            rows.append((cols[0], int(cols[1])))
    return rows


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    sets_path, tidk_dir = sys.argv[1], Path(sys.argv[2])
    outfile = sys.argv[3] if len(sys.argv) > 3 else None
    vp_outfile = sys.argv[4] if len(sys.argv) > 4 else None

    rows = list(csv.DictReader(open(sets_path), delimiter="\t"))
    identical = [r for r in rows if int(r["n_distinct"]) == 1 and 6 <= int(r["core_length"]) <= 15]

    out_rows = []
    for r in identical:
        sp, length, dominant = r["species"], int(r["core_length"]), r["dominant"]
        tidk_path = tidk_dir / f"{sp}.tidk.tsv"
        if not tidk_path.exists():
            continue
        candidates = load_tidk(tidk_path)
        dom_count = next((c for seq, c in candidates if seq == dominant), 0)

        best = None
        for seq, count in candidates:
            if seq == dominant or len(seq) != length or count < MIN_SECOND_COUNT:
                continue
            diff_n, shift = best_rotation(dominant, seq)
            if diff_n == 0:
                continue  # pure rotation of the same sequence, not a real second variant
            if diff_n > MAX_DIFF:
                continue
            if best is None or count > best[1]:
                best = (seq, count, diff_n, shift)

        if best is None:
            continue
        second, second_count, diff_n, shift = best
        out_rows.append(dict(
            species=sp, core_length=length, n_loci=r["n_loci"],
            dominant=dominant, dominant_tidk_count=dom_count,
            second=second, second_tidk_count=second_count,
            best_diff_n=diff_n, best_shift=shift,
        ))

    cols = ["species", "core_length", "n_loci", "dominant", "dominant_tidk_count",
            "second", "second_tidk_count", "best_diff_n", "best_shift"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in out_rows:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    print(f"[find_identical_set_second_candidates] {len(identical)} identical sets checked, "
          f"{len(out_rows)} have a viable close second tidk candidate (independent of TR)",
          file=sys.stderr)

    if vp_outfile:
        vp_cols = ["species", "core_length", "n_loci", "n_distinct", "variant", "n_copies",
                   "is_modal", "target_example", "raw_diff_n", "best_shift", "best_diff_n",
                   "is_pure_rotation", "diff_positions", "diff_bases"]
        with open(vp_outfile, "w") as vf:
            print("\t".join(vp_cols), file=vf)
            for r in out_rows:
                common = dict(species=r["species"], core_length=r["core_length"], n_loci=r["n_loci"],
                              n_distinct=2, target_example="", raw_diff_n=0, best_shift=0,
                              is_pure_rotation=0, diff_positions="", diff_bases="")
                dom = dict(common, variant=r["dominant"], n_copies=r["n_loci"], is_modal=1, best_diff_n=0)
                sec = dict(common, variant=r["second"], n_copies=1, is_modal=0, best_diff_n=r["best_diff_n"])
                print("\t".join(str(dom[c]) for c in vp_cols), file=vf)
                print("\t".join(str(sec[c]) for c in vp_cols), file=vf)
        print(f"[find_identical_set_second_candidates] variant-positions adapter written to {vp_outfile}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
