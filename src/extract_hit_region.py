#!/usr/bin/env python3
"""Extract a genomic region (e.g. a tblastn hit) and translate it in the
reading frame with the fewest stop codons, for use as a new curated TERT
query in inputs/plant_TERT_regions.fasta / .faa.

Intended for regions identified from outputs/tert_tblastn/<species>.tbl
that have no gene annotation to fall back on (e.g. Santalales), so the
frame has to be inferred rather than read off a GFF.

Usage:
  extract_hit_region.py <genome.fa[.gz]> <chrom> <start> <end> <strand> <name>

  start/end : 1-based inclusive, in forward-strand genome coordinates
              (i.e. min/max of BLAST sstart/send, regardless of hit strand)
  strand    : + or - (use - when the tblastn hit had sstart > send)
  name      : short identifier used in the FASTA headers

Prints two FASTA records to stdout: the nucleotide region, then the
translated protein in the best frame, with the stop-codon count of every
frame reported to stderr so the choice is auditable.
"""
import gzip
import sys

from Bio import SeqIO
from Bio.Seq import Seq


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    genome_path, chrom, start, end, strand, name = sys.argv[1:]
    start, end = int(start), int(end)

    with open_maybe_gz(genome_path) as fh:
        record = next(r for r in SeqIO.parse(fh, "fasta") if r.id == chrom)

    region = record.seq[start - 1:end]
    if strand == "-":
        region = region.reverse_complement()
    elif strand != "+":
        sys.exit(f"strand must be + or -, got {strand!r}")

    best_frame, best_prot, best_stops = None, None, None
    for frame in range(3):
        trimmed = region[frame:]
        trimmed = trimmed[: len(trimmed) - len(trimmed) % 3]
        prot = trimmed.translate()
        stops = prot.count("*")
        print(f"frame {frame}: {stops} stop codon(s)", file=sys.stderr)
        if best_stops is None or stops < best_stops:
            best_frame, best_prot, best_stops = frame, prot, stops

    print(
        f"[extract_hit_region] chose frame {best_frame} ({best_stops} stop codons)",
        file=sys.stderr,
    )

    header_common = (
        f"{chrom}:{start}-{end}({strand}) genomic hit region, frame {best_frame}"
    )
    print(f">{name} {header_common}")
    print(str(region))
    print(f">{name} {header_common} [translated]")
    print(str(best_prot))


if __name__ == "__main__":
    main()
