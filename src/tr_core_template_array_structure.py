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

Both questions are answered TWICE per set: once over the full extraction
window, and once restricted to "telomere-proper" (the contiguous run of
matches starting exactly at the true chromosome terminus, up to the
first large gap - see telomere_proper_prefix()'s docstring). A fixed
extraction window can include subtelomeric/interstitial sequence beyond
the true telomere, and interstitial telomeric repeats are known to be
more degenerate than terminal ones (Belyayev et al. 2023) - the
full_*/proper_* column pairs in the output let you see how much this
matters (in practice: a lot - several sets lose their entire minor-
variant signal once restricted).

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


import os
GAP_THRESHOLD = int(os.environ.get("GAP_THRESHOLD", 50))  # bp of unmatched sequence that
# ends the "telomere-proper" run from the terminus - overridable via env var for a
# sensitivity check (see notes/tr_core_template_array_structure.md's caveats)


def telomere_proper_prefix(track, end):
    """Walking from the true chromosome terminus (position 0 of the track
    for 5prime; the LAST tile for 3prime, since extract_terminal_prepared
    returns a suffix ending at the chromosome's actual last base), return
    the sub-list of track entries up to (not including) the first gap
    exceeding GAP_THRESHOLD bp. Plant telomere length varies enormously
    across species (documented range ~0.3kb to 200kb+ - see e.g. Watson &
    Riha 2010 FEBS Lett, Procházková Schrumpfová et al. 2016 Front Plant
    Sci), so a fixed 5000bp extraction window could be entirely within
    the telomere for a long-telomere species or well past it into
    subtelomeric/interstitial sequence for a short-telomere one - this
    function doesn't assume either: it just walks from the true terminus
    and stops at the first large break in matches, wherever that happens
    to be, rather than at a fixed distance. Interstitial telomeric
    repeats (ITRs) are known to be MORE degenerate than terminal arrays
    (Belyayev et al. 2023), so pooling in whatever lies beyond that break
    risks mistaking ITR noise for real terminal array structure. Caveat:
    GAP_THRESHOLD is a somewhat arbitrary tunable parameter - a real
    within-telomere transition between variant blocks could in principle
    exceed it by chance for a species with a long, genuinely degenerate
    telomere, wrongly truncating what's still bona fide telomere; not
    tested for sensitivity to this choice yet. Small gaps (<=GAP_THRESHOLD)
    are tolerated as normal degenerate units within an otherwise-dense array."""
    ordered = track if end == "5prime" else list(reversed(track))
    proper = []
    cumulative_gap = 0
    for entry in ordered:
        start, tend, label = entry
        if label is None:
            cumulative_gap += (tend - start)
            if cumulative_gap > GAP_THRESHOLD:
                break
            proper.append(entry)
        else:
            cumulative_gap = 0
            proper.append(entry)
    return proper


def block_structure(per_chrom_labels):
    max_block = 0
    best_chrom, best_stats = None, dict(n0=0, n1=0, n_runs=None, expected=None, z=None)
    for (chrom, end), labels in per_chrom_labels.items():
        max_block = max(max_block, longest_run_of_ones(labels))
        if labels.count(1) > best_stats["n1"]:
            best_stats = runs_test(labels)
            best_chrom = f"{chrom}:{end}"
    return max_block, best_chrom, best_stats


def verdict_for(n_minor, max_block, best_stats):
    if n_minor == 0:
        return "NO_VARIATION_OBSERVED"
    elif max_block >= 2:
        return f"CLUSTERED (block of {max_block} contiguous minor-variant units found)"
    elif best_stats["z"] is not None and best_stats["z"] < -1.96:
        return "CLUSTERED (runs test)"
    elif best_stats["z"] is not None and best_stats["z"] > 1.96:
        return "REGULAR/DISPERSED (runs test)"
    else:
        return "CONSISTENT WITH RANDOM (isolated single-unit substitutions)"


def analyze_species(species, species_sets, positional_dir, window):
    fa, tmpdir = prepare_genome(species)
    results = []
    try:
        for (sp, core_length), s in species_sets:
            modal, variants = s["modal"], s["variants"]
            targets = all_main_chromosomes(sp, modal, positional_dir)
            # raw per-chromosome tile sequence (base variant string per matched
            # tile, in genome order) - kept so the structural split can be
            # decided AFTER seeing which variant actually dominates the
            # telomere (not assumed to be the TR-copy modal - see LABEL_FLIP).
            # Computed twice: over the FULL extraction window, and restricted
            # to the terminus-contiguous "telomere-proper" prefix (see
            # telomere_proper_prefix's docstring) - a fixed 5000bp window is
            # much longer than a typical plant telomere, so full-window
            # structure can be dominated by subtelomeric/interstitial noise.
            full_counts, proper_counts = Counter(), Counter()
            per_chrom_full, per_chrom_proper = {}, {}
            for chrom, end, approx_length in targets:
                seq, offset, length = extract_terminal_prepared(fa, chrom, end, window, approx_length=approx_length)
                track = tile(seq, variants)
                full_labels = [strip_label(label) for _, _, label in track if label]
                proper_track = telomere_proper_prefix(track, end)
                proper_labels = [strip_label(label) for _, _, label in proper_track if label]
                for base in full_labels:
                    full_counts[base] += 1
                for base in proper_labels:
                    proper_counts[base] += 1
                if full_labels:
                    per_chrom_full[(chrom, end)] = full_labels
                if proper_labels:
                    per_chrom_proper[(chrom, end)] = proper_labels

            def summarize(counts, per_chrom):
                # structure test: telomere-majority variant vs everything else
                # pooled (not TR-copy modal vs rest - those disagree exactly in
                # the cases, like LABEL_FLIP, where the structure question is
                # most interesting, and using the wrong reference makes an
                # enormous majority block look like a spurious "cluster").
                dominant, dom_n = counts.most_common(1)[0] if counts else (None, 0)
                n_modal = counts.get(modal, 0)
                n_minor = sum(v for k, v in counts.items() if k != dominant)
                labels = {k: [0 if b == dominant else 1 for b in v] for k, v in per_chrom.items()}
                max_block, best_chrom, best_stats = block_structure(labels)
                verdict = verdict_for(n_minor, max_block, best_stats)
                return dict(
                    dominant=dominant, dom_n=dom_n, n_modal=n_modal, n_minor=n_minor,
                    copy_dominant_is_dominant=int(modal == dominant),
                    n_chroms_with_minor=sum(1 for lbl in labels.values() if 1 in lbl),
                    max_block=max_block, runs_chrom=best_chrom or "",
                    runs_n0=best_stats["n0"], runs_n1=best_stats["n1"],
                    runs_observed=best_stats["n_runs"], runs_expected=best_stats["expected"],
                    runs_z=best_stats["z"], verdict=verdict,
                )

            full = summarize(full_counts, per_chrom_full)
            proper = summarize(proper_counts, per_chrom_proper)

            results.append(dict(
                species=sp, core_length=core_length, n_loci=s["n_loci"],
                n_distinct_variants=len(variants), modal=modal,
                n_chroms_checked=len({c for c, _, _ in targets}),
                full_dominant=full["dominant"], full_dom_n=full["dom_n"],
                full_n_minor=full["n_minor"], full_copy_is_dominant=full["copy_dominant_is_dominant"],
                full_max_block=full["max_block"], full_runs_z=full["runs_z"], full_verdict=full["verdict"],
                proper_dominant=proper["dominant"], proper_dom_n=proper["dom_n"],
                proper_n_minor=proper["n_minor"], proper_copy_is_dominant=proper["copy_dominant_is_dominant"],
                proper_max_block=proper["max_block"], proper_runs_z=proper["runs_z"], proper_verdict=proper["verdict"],
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

    cols = ["species", "core_length", "n_loci", "n_distinct_variants", "modal", "n_chroms_checked",
            "full_dominant", "full_dom_n", "full_n_minor", "full_copy_is_dominant",
            "full_max_block", "full_runs_z", "full_verdict",
            "proper_dominant", "proper_dom_n", "proper_n_minor", "proper_copy_is_dominant",
            "proper_max_block", "proper_runs_z", "proper_verdict"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in all_results:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()

    n_flip = sum(1 for r in all_results if not r["proper_copy_is_dominant"])
    n_full_clustered = sum(1 for r in all_results if r["full_verdict"].startswith("CLUSTERED"))
    n_proper_clustered = sum(1 for r in all_results if r["proper_verdict"].startswith("CLUSTERED"))
    n_proper_absent = sum(1 for r in all_results if r["proper_verdict"] == "NO_VARIATION_OBSERVED")
    print(f"\n[tr_core_template_array_structure] {len(all_results)} sets analyzed; "
          f"{n_flip} have a telomere-proper-dominant variant different from the TR-copy-count "
          f"dominant; full-window clustered: {n_full_clustered}; telomere-proper-restricted "
          f"clustered: {n_proper_clustered}; telomere-proper shows no minor variant at all: "
          f"{n_proper_absent}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
