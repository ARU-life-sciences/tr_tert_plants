#!/usr/bin/env python3
"""Checks the "orientation ambiguity" flagged in tr_repeat_correlation.md:
85% of TR Template-to-genome-repeat matches are in the `direct` sense
rather than the textbook `revcomp` templating-mechanism relationship
(RNA template -> reverse-transcribed, complementary G-strand DNA).

Two checks:
1. Strand confound: does match_orientation (direct/revcomp) correlate
   with which genomic strand the TR gene itself sits on (+/-)? If it did,
   that would point to a strand-handling bug in extraction. It doesn't
   (checked below) - ruling out that explanation.
2. C-rich/G-rich composition of `direct` matches: if both the extracted
   Template and tidk's reported repeat string are independently,
   correctly C-rich (as expected - Fajkus et al. 2019's own
   experimentally-validated AtTR template, CTAAACCCT, is entirely C-rich,
   since the RNA template is naturally the complementary-strand sequence
   to the G-rich DNA it guides synthesis of), then a literal match
   between two C-rich strings is `direct` by construction - this would
   resolve the "ambiguity" as a strand-convention artifact (tidk
   reporting the lexicographically-smallest rotation of a
   strand-ambiguous repeat family, which for the canonical TTTAGGG/
   AAACCCT family happens to be the C-rich AAACCCT), not a genuine
   contradiction of the templating mechanism.

Usage:
  check_orientation_ambiguity.py <all_species.tsv>
"""
import csv
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def rotations(s):
    return [s[i:] + s[:i] for i in range(len(s))]


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))

    print("=== Check 1: match_orientation vs gene strand (confound test) ===")
    from collections import Counter
    combo = Counter((r["strand"], r["match_orientation"]) for r in rows)
    for k, v in sorted(combo.items()):
        print(f"  strand={k[0]} orientation={k[1]}: {v}")
    print("  -> if not confounded, direct/revcomp ratio should be similar for + and -\n")

    print("=== Check 2: base composition of `direct` matches ===")
    direct = [r for r in rows if r["match_orientation"] == "direct" and int(r["match_len"]) > 0]
    g_rich = c_rich = tied = 0
    for r in direct:
        chunk = r["match_chunk"]
        g, c = chunk.count("G"), chunk.count("C")
        if g > c:
            g_rich += 1
        elif c > g:
            c_rich += 1
        else:
            tied += 1
    print(f"  total direct matches: {len(direct)}")
    print(f"  C-rich (more C than G): {c_rich} ({100*c_rich/len(direct):.1f}%)")
    print(f"  G-rich (more G than C): {g_rich} ({100*g_rich/len(direct):.1f}%)")
    print(f"  tied: {tied}\n")

    print("=== Check 3: is the C-rich rotation the lexicographic minimum? ===")
    canonical_Grich = "TTTAGGG"
    crich = revcomp(canonical_Grich)
    all_rotations = rotations(canonical_Grich) + rotations(crich)
    smallest = sorted(all_rotations)[0]
    print(f"  G-rich canonical: {canonical_Grich}, C-rich (revcomp): {crich}")
    print(f"  lexicographically smallest of all 14 rotations (both strands): {smallest}")
    print(f"  matches Fajkus et al. 2019's validated AtTR template family "
          f"(CTAAACCCT, C-rich)? {'yes' if smallest in crich + crich else 'no'}")


if __name__ == "__main__":
    main()
