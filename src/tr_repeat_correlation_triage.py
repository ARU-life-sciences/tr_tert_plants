#!/usr/bin/env python3
"""Triage the per-locus output of tr_repeat_correlation.py into a
shortlist worth positional confirmation (confirm_repeat_positional.py).

A locus is a candidate if its best match is (a) statistically significant
against the shuffle-composition null (shuffle_p < P_THRESHOLD), and (b)
not just a re-discovery of the canonical TTTAGGG repeat (that's the
expected/uninteresting case - most loci do this). Candidates are then
grouped per species and scored by:
  - how many independent loci in that species share the exact same
    matched chunk (a shared, specific, non-generic sequence recurring at
    multiple paralogous loci is strong evidence it's a real repeat
    variant, not one-off extraction noise)
  - how many nucleotides different that chunk's best 7bp window is from
    the nearest canonical rotation (min_edit_to_canonical) - low values
    (1) are plausibly just a degenerate/wobble position within an
    otherwise-canonical-type telomere; higher values are more novel
  - whether the matched chunk is low-complexity (<=2 distinct bases -
    e.g. poly-A/T) and therefore a weaker candidate even if "significant",
    since low-complexity strings are more prone to spurious long matches
    in composition-matched shuffle tests than genuinely specific motifs

Usage:
  tr_repeat_correlation_triage.py <all_species.tsv> [outfile]

Writes a per-species summary table (species, n_significant_loci,
n_sharing_top_chunk, top_chunk, top_chunk_genome_count,
min_edit_to_canonical, low_complexity, recommended) to stdout or outfile.
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_tr_domains import find_repeat_in_window

P_THRESHOLD = 0.05
MIN_SHARED_FOR_RECOMMEND = 3
MAX_LOW_COMPLEXITY_DISTINCT_BASES = 2

_CANON_BASE = "CCCTAAA"
_CANON_ROTATIONS = [_CANON_BASE[i:] + _CANON_BASE[:i] for i in range(len(_CANON_BASE))]


def is_canonical_chunk(chunk):
    if len(chunk) < 6:
        return False
    return bool(find_repeat_in_window(chunk + chunk, "TTTAGGG")) or bool(find_repeat_in_window(chunk, "TTTAGGG"))


def is_low_complexity(chunk, max_distinct=MAX_LOW_COMPLEXITY_DISTINCT_BASES):
    return len(set(chunk)) <= max_distinct


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def min_edit_to_canonical(chunk):
    """Best (lowest) Hamming distance between any 7bp window of chunk
    (tandem-doubled, to catch wrap-around windows) and any canonical 7bp
    rotation."""
    seq = chunk + chunk
    best = 99
    for i in range(len(chunk)):
        window = seq[i:i + 7]
        if len(window) < 7:
            continue
        for rot in _CANON_ROTATIONS:
            best = min(best, hamming(window, rot))
    return best


def load_rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            yield dict(zip(header, cols))


def triage(path):
    species_chunks = defaultdict(list)
    for row in load_rows(path):
        mlen = int(row["match_len"])
        p = float(row["shuffle_p"])
        chunk = row["match_chunk"]
        if mlen == 0 or p >= P_THRESHOLD:
            continue
        if is_canonical_chunk(chunk):
            continue
        species_chunks[row["species"]].append(row)

    results = []
    for sp, entries in species_chunks.items():
        chunks = [e["match_chunk"] for e in entries]
        cc = Counter(chunks)
        top_chunk, shared = cc.most_common(1)[0]
        top_count = max(int(e["tidk_count"]) for e in entries if e["match_chunk"] == top_chunk)
        dist = min_edit_to_canonical(top_chunk)
        low_complexity = is_low_complexity(top_chunk)
        recommended = shared >= MIN_SHARED_FOR_RECOMMEND and not low_complexity
        results.append(dict(
            species=sp, n_significant_loci=len(entries), n_sharing_top_chunk=shared,
            top_chunk=top_chunk, top_chunk_genome_count=top_count,
            min_edit_to_canonical=dist, low_complexity=low_complexity,
            recommended_for_positional_confirmation=recommended,
        ))
    results.sort(key=lambda r: (r["min_edit_to_canonical"], -r["n_sharing_top_chunk"], -r["top_chunk_genome_count"]))
    return results


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else None
    results = triage(path)

    cols = ["species", "n_significant_loci", "n_sharing_top_chunk", "top_chunk",
            "top_chunk_genome_count", "min_edit_to_canonical", "low_complexity",
            "recommended_for_positional_confirmation"]
    out = open(outfile, "w") if outfile else sys.stdout
    print("\t".join(cols), file=out)
    for r in results:
        print("\t".join(str(r[c]) for c in cols), file=out)
    if outfile:
        out.close()


if __name__ == "__main__":
    main()
