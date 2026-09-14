#!/usr/bin/env python3
"""Control-group counterpart to tr_core_template_array_structure.py: for
species with FULLY IDENTICAL TR gene copies (no paralog sequence
divergence at all) and a plausible independent second repeat candidate
(from find_identical_set_second_candidates.py - discovered from that
species' own tidk-explore output, blind to TR gene content), exhaustively
walk every main chromosome, both ends, and test for the same
blocky/clustered/regular-alternation array structure found in the 44
TR-paralog-differing sets.

Differs from tr_core_template_array_structure.py only in how "main
chromosome" is defined: the differing-set pipeline reused an existing
tidk-search window-count threshold (confirm_repeat_positional.py's
output, already computed for other reasons). Here there is no such
pre-existing scan, and re-running tidk search per species just to
enumerate chromosomes is unnecessary compute for something the genome
assembly itself already answers directly - main chromosome = long enough
scaffold, judged from the genome's own sequence lengths (samtools faidx),
not from how much of any one candidate repeat happens to show up there.
This is arguably a more principled and consistent definition than the
window-count threshold used for the differing sets, not just a shortcut;
noted as a methods difference between the two groups, not assumed
equivalent.

Usage:
  tr_identical_set_array_structure.py <identical_set_variant_positions.tsv> [window_bp] [outfile]
"""
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from walk_terminal_repeat_array import prepare_genome, extract_terminal_prepared, tile
from tr_core_template_array_structure import (
    load_sets, strip_label, runs_test, longest_run_of_ones,
    telomere_proper_prefix, block_structure, verdict_for, GAP_THRESHOLD,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WINDOW = 5000
SAMTOOLS = "/software/team301/samtools/samtools"
MIN_MAIN_LENGTH = 500_000  # matches MIN_MAIN_WINDOWS=100 * TIDK_WINDOW_BP=5000 used for the differing sets
MAX_MAIN_CHROMS = 40       # matches the differing-set pipeline's cap


def main_chromosomes_from_fasta(fa_path):
    """Every chromosome at least MIN_MAIN_LENGTH bp, longest first, capped
    at MAX_MAIN_CHROMS - both ends each. Uses samtools faidx directly on
    the already-decompressed genome (no repeat-specific scan needed)."""
    subprocess.run([SAMTOOLS, "faidx", str(fa_path)], check=True, capture_output=True)
    fai = Path(str(fa_path) + ".fai")
    sized = []
    with open(fai) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            name, length = cols[0], int(cols[1])
            if length >= MIN_MAIN_LENGTH:
                sized.append((name, length))
    sized.sort(key=lambda kv: -kv[1])
    sized = sized[:MAX_MAIN_CHROMS]

    out = []
    for chrom, length in sized:
        out.append((chrom, "5prime", None))
        out.append((chrom, "3prime", length))
    return out


def analyze_species(species, species_sets, window):
    fa, tmpdir = prepare_genome(species)
    results = []
    try:
        targets = main_chromosomes_from_fasta(fa)
        for (sp, core_length), s in species_sets:
            modal, variants = s["modal"], s["variants"]
            from collections import Counter
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
    window = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_WINDOW
    outfile = sys.argv[3] if len(sys.argv) > 3 else None

    sets = load_sets(variant_positions_path)
    from collections import defaultdict
    by_species = defaultdict(list)
    for key, s in sets.items():
        by_species[key[0]].append((key, s))

    all_results = []
    for i, (species, species_sets) in enumerate(sorted(by_species.items()), 1):
        print(f"[tr_identical_set_array_structure] ({i}/{len(by_species)}) {species}: "
              f"{len(species_sets)} set(s)", file=sys.stderr)
        try:
            all_results.extend(analyze_species(species, species_sets, window))
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)

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

    n_full_clustered = sum(1 for r in all_results if r["full_verdict"].startswith("CLUSTERED"))
    n_proper_clustered = sum(1 for r in all_results if r["proper_verdict"].startswith("CLUSTERED"))
    n_proper_absent = sum(1 for r in all_results if r["proper_verdict"] == "NO_VARIATION_OBSERVED")
    print(f"\n[tr_identical_set_array_structure] {len(all_results)} sets analyzed; "
          f"full-window clustered: {n_full_clustered}; telomere-proper-restricted clustered: "
          f"{n_proper_clustered}; telomere-proper shows no minor variant at all: {n_proper_absent}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
