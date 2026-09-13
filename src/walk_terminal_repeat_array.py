#!/usr/bin/env python3
"""Base-pair-resolved answer to a question tidk can't answer: given two (or
more) candidate repeat-unit variants, walk the ACTUAL terminal genomic
sequence of a chromosome and report, position by position, which variant
occurs where - not just how many times each occurs somewhere in a window.

tidk (explore or search) counts occurrences of a literal string per
window; it cannot tell you whether two variants are interspersed within
one repeat array, segregated into separate blocks, or arranged in any
other specific order. This does, by greedily tiling the raw sequence
left-to-right and reporting which variant (if any) matched at each
position, and rendering the result as a linear track.

Usage:
  walk_terminal_repeat_array.py <species> <chromosome> <end> <variant1> [variant2 ...] [--window BP] [--outdir DIR]

  end       5prime or 3prime - which chromosome end to pull sequence from
  variantN  exact repeat-unit strings to scan for (any length); reverse
            complement is also tried automatically
  --window  bp of terminal sequence to pull (default 3000)
"""
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ESL_SFETCH = "/software/team301/hmmer-3.4/easel/miniapps/esl-sfetch"
DEFAULT_WINDOW = 3000

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def load_genome_path(species):
    with open(PROJECT_ROOT / "inputs" / "dtol_plant_paths.txt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == species:
                return parts[1]
    return None


def chrom_length(fa_path, chrom):
    out = subprocess.run([ESL_SFETCH, "--index", str(fa_path)], capture_output=True)
    fai = Path(str(fa_path) + ".ssi")
    # esl-sfetch has no direct "length" query; fetch whole seq via -n dummy trick isn't ideal,
    # so instead use a large end coordinate and let esl-sfetch clamp - simplest: fetch with -c 1..0 (whole seq)
    result = subprocess.run([ESL_SFETCH, "-c", "1..0", str(fa_path), chrom],
                             check=True, capture_output=True, text=True)
    seq = "".join(result.stdout.splitlines()[1:])
    return len(seq), seq


def extract_terminal(species, chrom, end, window):
    genome_path = load_genome_path(species)
    if genome_path is None:
        sys.exit(f"species {species!r} not in dtol_plant_paths.txt")
    tmpdir = Path(tempfile.mkdtemp(prefix=f"walk_{species}_"))
    try:
        fa = tmpdir / "genome.fa"
        if genome_path.endswith(".gz"):
            with gzip.open(genome_path, "rt") as fin, open(fa, "w") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copyfile(genome_path, fa)
        subprocess.run([ESL_SFETCH, "--index", str(fa)], check=True, capture_output=True)
        length, full_seq = chrom_length(fa, chrom)
        if end == "5prime":
            seq = full_seq[:window]
            offset = 0
        elif end == "3prime":
            seq = full_seq[-window:]
            offset = length - len(seq)
        else:
            sys.exit("end must be 5prime or 3prime")
        return seq, offset, length
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def tile(seq, variants):
    """Greedy left-to-right tiling: at each untiled position, try every
    variant (and its revcomp) longest-first, take the first exact match,
    advance past it; if nothing matches, advance by 1 and mark as gap."""
    variant_forms = []
    for v in variants:
        variant_forms.append((v, v))
        rc = revcomp(v)
        if rc != v:
            variant_forms.append((rc, "revcomp(" + v + ")"))
    variant_forms.sort(key=lambda x: -len(x[0]))

    i = 0
    n = len(seq)
    track = []  # (start, end, label) where label is variant string or None for gap
    while i < n:
        matched = None
        for form, label in variant_forms:
            if seq[i:i + len(form)] == form:
                matched = (i, i + len(form), label)
                break
        if matched:
            track.append(matched)
            i = matched[1]
        else:
            track.append((i, i + 1, None))
            i += 1
    return track


def render(track, seq, offset):
    """Compact run-length view: consecutive tiles with the same label
    collapsed to one entry with a repeat count and genomic position."""
    runs = []
    for start, end, label in track:
        if runs and runs[-1][2] == label:
            runs[-1] = (runs[-1][0], end, label, runs[-1][3] + 1)
        else:
            runs.append((start, end, label, 1))
    lines = []
    for start, end, label, n in runs:
        if label is None:
            if n <= 3:
                continue  # don't clutter with tiny 1-3bp gaps
            lines.append(f"  [gap {n}bp] genome_pos={offset+start+1}-{offset+end}")
        else:
            lines.append(f"  {label:20s} x{n:<4d} genome_pos={offset+start+1}-{offset+end}")
    return lines


def main():
    args = sys.argv[1:]
    if len(args) < 4:
        sys.exit(__doc__)
    species, chrom, end = args[0], args[1], args[2]
    rest = args[3:]
    window = DEFAULT_WINDOW
    variants = []
    i = 0
    while i < len(rest):
        if rest[i] == "--window":
            window = int(rest[i + 1]); i += 2
        else:
            variants.append(rest[i]); i += 1

    print(f"[walk_terminal_repeat_array] {species} {chrom} {end} window={window}bp variants={variants}", file=sys.stderr)
    seq, offset, length = extract_terminal(species, chrom, end, window)
    print(f"[walk_terminal_repeat_array] chromosome length={length}, extracted {len(seq)}bp at offset {offset}", file=sys.stderr)

    track = tile(seq, variants)
    counts = {}
    for _, _, label in track:
        if label:
            counts[label] = counts.get(label, 0) + 1
    print(f"\n=== {species} {chrom} {end} ({len(seq)}bp scanned) ===")
    print("Tile counts:", counts)
    print("\nLinear order (runs, 5'->3' of the extracted window):")
    for line in render(track, seq, offset):
        print(line)


if __name__ == "__main__":
    main()
