#!/usr/bin/env python3
"""Generalization of the Carlina vulgaris analysis
(notes/carlina_vulgaris_telomere_repeat.md,
notes/tr_repeat_correlation.md): for every TR-homologous locus in a
species, does its specific Template sequence contain a long, exact,
statistically-unlikely-by-chance substring match to one of that species'
own genome-wide (chromosome-terminal-proximal, per tidk explore's
--distance restriction) repeat units? If so, which one, and how
significant is the match against a composition-matched random-shuffle
null?

This tests a finer question than "is the Template canonical or not": does
the SPECIFIC sequence of a given TR paralog's Template correspond to a
SPECIFIC repeat variant actually present at this genome's chromosome
termini - i.e. does point-mutation-level variation between TR gene copies
track point-mutation-level variation in the telomeric repeat itself.

Inputs (already produced by the rest of the pipeline, no genome access
needed - this is pure text-file analysis, hence no bsub variant):
  outputs/tr_template_loci/<species>.tsv   (src/extract_tr_template_all_loci.py)
  outputs/tidk/<species>.tidk.tsv          (src/run_tidk.bash)

Usage:
  tr_repeat_correlation.py <species> [n_shuffles]     # single species, prints TSV to stdout
  tr_repeat_correlation.py --all [n_shuffles] [outfile]  # every multi-locus species

Output columns (TSV):
  species  target  coords  strand  evalue  template  match_orientation
  match_len  match_chunk  tidk_repeat  tidk_count  shuffle_p  n_shuffles
"""
import glob
import random
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
MIN_MATCH_LEN = 6
DEFAULT_N_SHUFFLES = 200

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def load_loci(species):
    """[(target, coords, strand, evalue, template), ...] for OK rows."""
    path = OUTPUT_DIR / "tr_template_loci" / f"{species}.tsv"
    loci = []
    if not path.exists():
        return loci
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 7 and cols[1] == "OK":
                loci.append((cols[2], cols[3], cols[4], cols[5], cols[6]))
    return loci


def load_tidk_repeats(species):
    """[(repeat_string, count), ...] full list, ranked as tidk reported."""
    path = OUTPUT_DIR / "tidk" / f"{species}.tidk.tsv"
    repeats = []
    if not path.exists():
        return repeats
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[1].isdigit():
                repeats.append((parts[0], int(parts[1])))
    return repeats


def best_match_for_query(query, tidk_repeats, min_len=MIN_MATCH_LEN):
    """Longest exact substring shared between `query` and any tidk repeat
    string (checked at every rotation/offset within each repeat, since
    tidk repeats are themselves already-tandemised strings e.g.
    'AAACTGAAACTGAAACTG...'). Ties on length broken by higher tidk count
    (prefer the more genome-abundant/terminal-enriched family). Returns
    (length, chunk, repeat, count) or None. Searches longest-first with
    early exit, since the common case (a real hit) resolves at a long
    length quickly and the worst case (no hit) is comparatively rare."""
    max_len = min(30, len(query))
    for length in range(max_len, min_len - 1, -1):
        best_here = None
        for repeat, count in tidk_repeats:
            if len(repeat) < length:
                continue
            for start in range(len(repeat) - length + 1):
                chunk = repeat[start:start + length]
                if chunk in query:
                    if best_here is None or count > best_here[2]:
                        best_here = (chunk, repeat, count)
        if best_here is not None:
            return (length, best_here[0], best_here[1], best_here[2])
    return None


def best_match_either_orientation(template, tidk_repeats, min_len=MIN_MATCH_LEN):
    """Try both the Template sequence as-is ('direct') and its reverse
    complement ('revcomp' - the orientation expected if this locus is
    acting as a conventional reverse-transcription template), return
    whichever gives the longer match, tagged with which orientation and
    the exact query string that matched (so the same query can be reused
    for the shuffle control)."""
    direct = best_match_for_query(template, tidk_repeats, min_len)
    rc_query = revcomp(template)
    rc = best_match_for_query(rc_query, tidk_repeats, min_len)
    candidates = []
    if direct:
        candidates.append(("direct", template, direct))
    if rc:
        candidates.append(("revcomp", rc_query, rc))
    if not candidates:
        return None
    candidates.sort(key=lambda c: (c[2][0], c[2][3]), reverse=True)
    return candidates[0]


def shuffle_pvalue(query, tidk_repeats, real_len, n_shuffles, rng, min_len=MIN_MATCH_LEN):
    """Empirical p-value: fraction of random same-composition shuffles of
    `query` whose best match to the same tidk repeat list is >= real_len."""
    chars = list(query)
    n_ge = 0
    for _ in range(n_shuffles):
        rng.shuffle(chars)
        shuffled = "".join(chars)
        m = best_match_for_query(shuffled, tidk_repeats, min_len)
        if m is not None and m[0] >= real_len:
            n_ge += 1
    return n_ge / n_shuffles


def analyse_species(species, n_shuffles=DEFAULT_N_SHUFFLES, seed=42):
    """Yields one result tuple per locus (only species with >=2 loci and
    tidk data produce anything)."""
    loci = load_loci(species)
    if len(loci) < 2:
        return
    tidk_repeats = load_tidk_repeats(species)
    if not tidk_repeats:
        return
    rng = random.Random(seed)
    for target, coords, strand, evalue, template in loci:
        best = best_match_either_orientation(template, tidk_repeats)
        if best is None:
            yield (species, target, coords, strand, evalue, template,
                   "none", 0, "", "", 0, 1.0, n_shuffles)
            continue
        orientation, query, (length, chunk, repeat, count) = best
        p = shuffle_pvalue(query, tidk_repeats, length, n_shuffles, rng)
        yield (species, target, coords, strand, evalue, template,
               orientation, length, chunk, repeat, count, p, n_shuffles)


def multi_locus_species(assembly_list=None):
    if assembly_list is None:
        assembly_list = PROJECT_ROOT / "inputs" / "dtol_plant_paths.txt"
    species_in_list = []
    with open(assembly_list) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                species_in_list.append(parts[0])
    in_list = set(species_in_list)
    out = []
    for f in sorted(glob.glob(str(OUTPUT_DIR / "tr_template_loci" / "*.tsv"))):
        sp = Path(f).stem
        if sp not in in_list:
            continue
        if len(load_loci(sp)) >= 2:
            out.append(sp)
    return out


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)

    header = ("species\ttarget\tcoords\tstrand\tevalue\ttemplate\t"
              "match_orientation\tmatch_len\tmatch_chunk\ttidk_repeat\t"
              "tidk_count\tshuffle_p\tn_shuffles")

    if args[0] == "--all":
        n_shuffles = int(args[1]) if len(args) > 1 else DEFAULT_N_SHUFFLES
        outfile = args[2] if len(args) > 2 else None
        species_list = multi_locus_species()
        out = open(outfile, "w") if outfile else sys.stdout
        print(header, file=out)
        for i, sp in enumerate(species_list, 1):
            print(f"[{i}/{len(species_list)}] {sp}", file=sys.stderr)
            for row in analyse_species(sp, n_shuffles=n_shuffles):
                print("\t".join(str(x) for x in row), file=out)
            if outfile:
                out.flush()
        if outfile:
            out.close()
    else:
        species = args[0]
        n_shuffles = int(args[1]) if len(args) > 1 else DEFAULT_N_SHUFFLES
        print(header)
        for row in analyse_species(species, n_shuffles=n_shuffles):
            print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
