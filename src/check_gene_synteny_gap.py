#!/usr/bin/env python3
"""Checks whether a gene apparently missing from one species' genome
(by HMM/tblastn search - no hit) is a genuine sequence absence or an
assembly-completeness artifact ("assembly gap"), by testing synteny
against a close relative that DOES have a confident hit for the gene.

Extracts a window around the gene in the relative (reference) species,
blastn's it against the full genome of the species missing the gene, and
reports: (1) whether a confident, single-locus syntenic anchor exists
nearby (not just scattered repetitive-element noise), (2) whether the
gene body itself has any detectable hit at all, and (3) whether the
anchor sits near a scaffold end or an N-gap in the query-missing
species' assembly (which would reopen "assembly gap" as the explanation)
or safely mid-scaffold in fully resolved sequence (which rules it out).

Usage:
  check_gene_synteny_gap.py <reference_species> <ref_chrom> <ref_start> <ref_end> <missing_species> [outdir]

  ref_start/ref_end: 1-based coordinates in the reference species'
  genome, spanning a window comfortably wider than the gene itself
  (enough flanking sequence on both sides to find a real anchor).

Writes <outdir>/<reference>_<missing>_synteny_blastn.tsv (raw blastn
hits, outfmt 6 + evalue/bitscore) and prints a summary to stderr:
merged query-coverage intervals, gap regions with no hits, and an
N-content / scaffold-position check around the nearest anchor.
"""
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from walk_terminal_repeat_array import prepare_genome

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ESL_SFETCH = "/software/team301/hmmer-3.4/easel/miniapps/esl-sfetch"
SAMTOOLS = "/software/team301/samtools/samtools"
BLASTN = "/software/team301/ncbi-blast-2.16.0+/bin/blastn"
MAKEBLASTDB = "/software/team301/ncbi-blast-2.16.0+/bin/makeblastdb"
EVALUE = "1e-10"


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    ref_species, ref_chrom, ref_start, ref_end, missing_species = sys.argv[1:6]
    outdir = Path(sys.argv[6]) if len(sys.argv) > 6 else PROJECT_ROOT / "outputs" / "synteny_checks"
    outdir.mkdir(parents=True, exist_ok=True)

    ref_fa, ref_tmp = prepare_genome(ref_species)
    missing_fa, missing_tmp = prepare_genome(missing_species)
    try:
        query_fa = outdir / f"{ref_species}_{ref_chrom}_{ref_start}-{ref_end}.fa"
        result = subprocess.run(
            [ESL_SFETCH, "-c", f"{ref_start}..{ref_end}", str(ref_fa), ref_chrom],
            check=True, capture_output=True, text=True)
        query_fa.write_text(result.stdout)
        print(f"[check_gene_synteny_gap] extracted {ref_species} {ref_chrom}:{ref_start}-{ref_end} "
              f"-> {query_fa}", file=sys.stderr)

        db_prefix = outdir / f"{missing_species}_db"
        subprocess.run([MAKEBLASTDB, "-in", str(missing_fa), "-dbtype", "nucl",
                         "-out", str(db_prefix)], check=True, capture_output=True)

        hits_tsv = outdir / f"{ref_species}_{missing_species}_synteny_blastn.tsv"
        subprocess.run(
            [BLASTN, "-query", str(query_fa), "-db", str(db_prefix),
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
             "-evalue", EVALUE, "-out", str(hits_tsv)],
            check=True, capture_output=True)
        print(f"[check_gene_synteny_gap] blastn hits -> {hits_tsv}", file=sys.stderr)

        rows = []
        for line in hits_tsv.read_text().splitlines():
            c = line.split("\t")
            qs, qe = int(c[6]), int(c[7])
            rows.append((min(qs, qe), max(qs, qe), c[1], float(c[2]), int(c[8]), int(c[9]), float(c[10])))
        rows.sort()

        merged = []
        for s, e, *_ in rows:
            if merged and s <= merged[-1][1] + 50:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))

        print(f"\n[check_gene_synteny_gap] {len(rows)} raw hits, merged into "
              f"{len(merged)} covered intervals (query-relative coords):", file=sys.stderr)
        for s, e in merged:
            ref_coord = int(ref_start) + s
            print(f"  query {s}-{e} (={ref_species} {ref_chrom}:{ref_coord}..{int(ref_start)+e}), "
                  f"len={e-s}", file=sys.stderr)

        window_len = int(ref_end) - int(ref_start)
        covered = set()
        for s, e in merged:
            covered.update(range(s, e))
        gaps = []
        in_gap = False
        gap_start = None
        for i in range(window_len):
            if i not in covered and not in_gap:
                in_gap, gap_start = True, i
            elif i in covered and in_gap:
                gaps.append((gap_start, i))
                in_gap = False
        if in_gap:
            gaps.append((gap_start, window_len))
        print(f"\n[check_gene_synteny_gap] uncovered (no-hit) query intervals >=500bp:", file=sys.stderr)
        for s, e in gaps:
            if e - s >= 500:
                print(f"  query {s}-{e} (len={e-s}) - no blastn hit anywhere in {missing_species}",
                      file=sys.stderr)

        # scaffold-position / N-content check around the single best (highest-identity) anchor
        if rows:
            best = max(rows, key=lambda r: r[3])
            s, e, sub, pid, ss, se, ev = best
            subprocess.run([SAMTOOLS, "faidx", str(missing_fa)], check=True, capture_output=True)
            fai_line = next((l for l in Path(str(missing_fa) + ".fai").read_text().splitlines()
                              if l.startswith(sub + "\t")), None)
            scaf_len = int(fai_line.split("\t")[1]) if fai_line else None
            anchor_pos = min(ss, se)
            print(f"\n[check_gene_synteny_gap] best anchor: query {s}-{e} vs {missing_species} "
                  f"{sub}:{ss}-{se} (pid={pid:.1f}%, evalue={ev:.2e})", file=sys.stderr)
            if scaf_len:
                frac = anchor_pos / scaf_len
                print(f"  {sub} length={scaf_len}, anchor at {frac*100:.1f}% through scaffold "
                      f"({'near an end - assembly gap plausible' if frac < 0.02 or frac > 0.98 else 'mid-scaffold - assembly gap unlikely'})",
                      file=sys.stderr)
                check_start = max(1, anchor_pos - 15000)
                check_end = min(scaf_len, anchor_pos + 15000)
                region = subprocess.run(
                    [ESL_SFETCH, "-c", f"{check_start}..{check_end}", str(missing_fa), sub],
                    check=True, capture_output=True, text=True).stdout
                seq = "".join(l.strip() for l in region.splitlines()[1:])
                n_frac = seq.upper().count("N") / len(seq) if seq else 1.0
                print(f"  N-content in {check_start}-{check_end} window around anchor: "
                      f"{100*n_frac:.2f}% ({'assembly gap present' if n_frac > 0.05 else 'fully resolved sequence'})",
                      file=sys.stderr)
    finally:
        shutil.rmtree(ref_tmp, ignore_errors=True)
        shutil.rmtree(missing_tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
