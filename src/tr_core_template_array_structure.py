#!/usr/bin/env python3
"""For every differing core-template set (33 species, 44 sets), exhaustively
walk every main chromosome's both ends (base-pair resolution, not tidk's
window counting - same machinery as resolve_weak_sensitivity.py) and
answer two things raw window-count data can't:

1. Which variant is actually dominant IN THE TELOMERE by direct sequence
   count - using ALL distinct core variants for the set (not just the
   top-2 dominant/minority tracked by core_template_sets.tsv), so this
   also catches cases where a third/fourth minor TR-copy variant turns
   out to be the real telomeric majority.
2. Is the non-modal (TR-copy-minority) sequence scattered as isolated
   single-unit substitutions through an otherwise-uniform array, or does
   it form contiguous multi-unit BLOCKS (candidate higher-order-repeat
   structure)? Answered two ways: directly (longest contiguous run of
   non-modal tiles found anywhere) and statistically (a Wald-Wolfowitz
   runs test against the null of random interspersion, run on the single
   chromosome-end with the most non-modal hits, since that's the best-
   powered available series).

One genome decompression per species, reused across every set and every
chromosome-end for that species.

Usage:
  tr_core_template_array_structure.py <core_template_variant_positions.tsv> [positional_dir] [window_bp] [outfile]
"""
import csv
import math
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from walk_terminal_repeat_array import prepare_genome, extract_terminal_prepared, tile
from resolve_weak_sensitivity import all_main_chromosomes

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WINDOW = 5000


def load_sets(path):
    """(species, core_length) -> dict(modal=str, variants=[str,...], n_loci=int)"""
    sets = defaultdict(lambda: dict(modal=None, variants=[], n_loci=0))
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["species"], row["core_length"])
            s = sets[key]
            s["n_loci"] = int(row["n_loci"])
            s["variants"].append(row["variant"])
            if row["is_modal"] == "1":
                s["modal"] = row["variant"]
    return sets


def strip_label(label):
    if label.startswith("revcomp(") and label.endswith(")"):
        return label[len("revcomp("):-1]
    return label


def runs_test(labels):
    """labels: list of 0/1. Wald-Wolfowitz runs test vs random interspersion."""
    n0, n1 = labels.count(0), labels.count(1)
    n = n0 + n1
    if n0 == 0 or n1 == 0 or n < 2:
        return dict(n0=n0, n1=n1, n_runs=None, expected=None, z=None)
    n_runs = 1 + sum(1 for i in range(1, n) if labels[i] != labels[i - 1])
    mean = 1 + (2 * n0 * n1) / n
    var_num = 2 * n0 * n1 * (2 * n0 * n1 - n)
    var_den = n * n * (n - 1)
    var = var_num / var_den if var_den else 0
    sd = math.sqrt(var) if var > 0 else 0
    z = (n_runs - mean) / sd if sd > 0 else None
    return dict(n0=n0, n1=n1, n_runs=n_runs, expected=round(mean, 2),
                z=round(z, 2) if z is not None else None)


def longest_run_of_ones(labels):
    best = cur = 0
    for v in labels:
        cur = cur + 1 if v == 1 else 0
        best = max(best, cur)
    return best


def analyze_species(species, species_sets, positional_dir, window):
    fa, tmpdir = prepare_genome(species)
    results = []
    try:
        for (sp, core_length), s in species_sets:
            modal, variants = s["modal"], s["variants"]
            targets = all_main_chromosomes(sp, modal, positional_dir)
            telomere_counts = Counter()
            # raw per-chromosome tile sequence (base variant string per matched
            # tile, in genome order) - kept so the structural split can be
            # decided AFTER seeing which variant actually dominates the
            # telomere (not assumed to be the TR-copy modal - see LABEL_FLIP).
            per_chrom_seq = {}
            for chrom, end, approx_length in targets:
                seq, offset, length = extract_terminal_prepared(fa, chrom, end, window, approx_length=approx_length)
                track = tile(seq, variants)
                seq_labels = [strip_label(label) for _, _, label in track if label]
                for base in seq_labels:
                    telomere_counts[base] += 1
                if seq_labels:
                    per_chrom_seq[(chrom, end)] = seq_labels

            telomere_dominant, telomere_dom_n = (
                telomere_counts.most_common(1)[0] if telomere_counts else (None, 0)
            )
            n_modal = telomere_counts.get(modal, 0)
            n_nonmodal = sum(v for k, v in telomere_counts.items() if k != modal)
            n_telomere_minor = sum(v for k, v in telomere_counts.items() if k != telomere_dominant)

            # structure test: telomere-majority variant vs everything else
            # pooled (not TR-copy modal vs rest - those disagree exactly in
            # the cases, like LABEL_FLIP, where the structure question is
            # most interesting, and using the wrong reference makes an
            # enormous majority block look like a spurious "cluster").
            per_chrom_labels = {
                key: [0 if base == telomere_dominant else 1 for base in seq_labels]
                for key, seq_labels in per_chrom_seq.items()
            }
            max_block = 0
            best_chrom, best_stats = None, dict(n0=0, n1=0, n_runs=None, expected=None, z=None)
            for (chrom, end), labels in per_chrom_labels.items():
                max_block = max(max_block, longest_run_of_ones(labels))
                if labels.count(1) > best_stats["n1"]:
                    best_stats = runs_test(labels)
                    best_chrom = f"{chrom}:{end}"

            if n_telomere_minor == 0:
                verdict = "NO_VARIATION_OBSERVED"
            elif max_block >= 2:
                verdict = f"CLUSTERED (block of {max_block} contiguous minor-variant units found)"
            elif best_stats["z"] is not None and best_stats["z"] < -1.96:
                verdict = "CLUSTERED (runs test)"
            elif best_stats["z"] is not None and best_stats["z"] > 1.96:
                verdict = "REGULAR/DISPERSED (runs test)"
            else:
                verdict = "CONSISTENT WITH RANDOM (isolated single-unit substitutions)"

            results.append(dict(
                species=sp, core_length=core_length, n_loci=s["n_loci"],
                n_distinct_variants=len(variants), modal=modal,
                telomere_dominant=telomere_dominant, telomere_dom_n=telomere_dom_n,
                n_modal=n_modal, n_nonmodal=n_nonmodal,
                copy_dominant_is_telomere_dominant=int(modal == telomere_dominant),
                n_chroms_checked=len({c for c, _, _ in targets}),
                n_chroms_with_minor=sum(1 for lbls in per_chrom_labels.values() if 1 in lbls),
                max_minor_block=max_block,
                runs_test_chrom=best_chrom or "",
                runs_n0=best_stats["n0"], runs_n1=best_stats["n1"],
                runs_observed=best_stats["n_runs"], runs_expected=best_stats["expected"],
                runs_z=best_stats["z"],
                verdict=verdict,
            ))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return results


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    variant_positions_path = sys.argv[1]
    positional_dir = sys.argv[2] if len(sys.argv) > 2 else str(PROJECT_ROOT / "outputs" / "tr_repeat_correlation" / "positional")
    window = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_WINDOW
    outfile = sys.argv[4] if len(sys.argv) > 4 else None

    sets = load_sets(variant_positions_path)
    by_species = defaultdict(list)
    for key, s in sets.items():
        by_species[key[0]].append((key, s))

    all_results = []
    for i, (species, species_sets) in enumerate(sorted(by_species.items()), 1):
        print(f"[tr_core_template_array_structure] ({i}/{len(by_species)}) {species}: "
              f"{len(species_sets)} set(s)", file=sys.stderr)
        all_results.extend(analyze_species(species, species_sets, positional_dir, window))

    cols = ["species", "core_length", "n_loci", "n_distinct_variants", "modal",
            "telomere_dominant", "telomere_dom_n", "n_modal", "n_nonmodal",
            "copy_dominant_is_telomere_dominant", "n_chroms_checked", "n_chroms_with_minor",
            "max_minor_block", "runs_test_chrom", "runs_n0", "runs_n1",
            "runs_observed", "runs_expected", "runs_z", "verdict"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in all_results:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_flip = sum(1 for r in all_results if not r["copy_dominant_is_telomere_dominant"])
    n_clustered = sum(1 for r in all_results if r["verdict"].startswith("CLUSTERED"))
    print(f"\n[tr_core_template_array_structure] {len(all_results)} sets analyzed; "
          f"{n_flip} have a telomere-dominant variant different from the TR-copy-count "
          f"dominant; {n_clustered} show clustered (block/HOR-like) non-modal structure",
          file=sys.stderr)


if __name__ == "__main__":
    main()
