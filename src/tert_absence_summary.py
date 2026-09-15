#!/usr/bin/env python3
"""Summarizes, for every species in the assembly list, whether TERT was
found by either detection method (HMM/nhmmer, protein BLAST/tblastn), to
surface candidate genuine TERT absences dataset-wide (not just the two
cases - Viscum_album, Epilobium_hirsutum - found so far by hand).

A species counts as having a TERT HMM hit if any of its 3 per-species
nhmmer .tbl files (isoformx1, catalytic_subunit, tert1) has at least one
non-comment line. A species counts as having a tblastn hit if its .tbl
file is non-empty (0-byte file = a completed, valid, genuine zero-hit
result - see worker_tblastn_tert.bash's atomic write, also used by
pipeline_status.bash).

"Missing TERT" here means: no hit by EITHER method. Missing by only one
method is reported separately since it's a much weaker signal (the
other, more sensitive method still found something).

Usage:
  tert_absence_summary.py [assembly_list] [outfile]

Output columns: species, hmm_hits, tblastn_hits, status
  status one of: MISSING_BOTH, MISSING_HMM_ONLY, MISSING_TBLASTN_ONLY, PRESENT
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def count_hmm_hits(species):
    spp_dir = PROJECT_ROOT / "outputs" / "tert_tbls" / species
    total = 0
    files = ["isoformx1", "catalytic_subunit", "tert1"]
    found = 0
    for name in files:
        f = spp_dir / f"{species}_{name}.tbl"
        if not f.exists():
            continue
        found += 1
        with open(f) as fh:
            for line in fh:
                if not line.startswith("#") and line.strip():
                    total += 1
    return total, found


def count_tblastn_hits(species):
    f = PROJECT_ROOT / "outputs" / "tert_tblastn" / f"{species}.tbl"
    if not f.exists():
        return None
    with open(f) as fh:
        return sum(1 for _ in fh)


def main():
    assembly_list = Path(sys.argv[1]) if len(sys.argv) > 1 else PROJECT_ROOT / "inputs" / "dtol_plant_paths.txt"
    outfile = sys.argv[2] if len(sys.argv) > 2 else None

    rows = []
    with open(assembly_list) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            species = line.split("\t")[0]

            hmm_hits, hmm_files_found = count_hmm_hits(species)
            tblastn_hits = count_tblastn_hits(species)

            hmm_missing = hmm_files_found == 3 and hmm_hits == 0
            tblastn_missing = tblastn_hits is not None and tblastn_hits == 0

            if hmm_files_found < 3 or tblastn_hits is None:
                status = "INCOMPLETE_DATA"
            elif hmm_missing and tblastn_missing:
                status = "MISSING_BOTH"
            elif hmm_missing:
                status = "MISSING_HMM_ONLY"
            elif tblastn_missing:
                status = "MISSING_TBLASTN_ONLY"
            else:
                status = "PRESENT"

            rows.append((species, hmm_hits, tblastn_hits if tblastn_hits is not None else "NA", status))

    out = open(outfile, "w") if outfile else sys.stdout
    print("species\thmm_hits\ttblastn_hits\tstatus", file=out)
    for r in rows:
        print("\t".join(str(x) for x in r), file=out)
    if outfile:
        out.close()

    from collections import Counter
    counts = Counter(r[3] for r in rows)
    print(f"\n[tert_absence_summary] {len(rows)} species:", file=sys.stderr)
    for k, v in counts.most_common():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
