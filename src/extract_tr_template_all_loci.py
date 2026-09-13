#!/usr/bin/env python3
"""For one species, extract the Template domain from EVERY decent TR-homologous
locus (not just the best hit, unlike extract_tr_domains.py), so within-species
variation across paralogous/duplicated loci can be assessed.

Unlike extract_tr_domains.py's linear hmmfrom/hmmto -> alifrom/alito
interpolation (fine for a single best hit, but interpolation noise could be
mistaken for real variation when comparing multiple loci against each
other), this extracts each full hit region and hmmaligns all loci for the
species together against the production HMM, then slices the Template
columns from that real alignment - the same method used to validate the
Carlina vulgaris finding (see notes/carlina_vulgaris_telomere_repeat.md).

Usage:
  extract_tr_template_all_loci.py <species> <genome_fa_gz> <nhmmer_tbl> <hmm> <outdir>

Writes <outdir>/<species>.tsv: one row per locus (target, coords, strand,
E-value, ungapped Template sequence), plus a summary line.
"""
import gzip
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ESL_SFETCH = "/software/team301/hmmer-3.4/easel/miniapps/esl-sfetch"
HMMALIGN = "/software/team301/hmmer-3.4/src/hmmalign"
EVALUE_THRESHOLD = 1e-4
MIN_SPAN = 100  # min hmm columns covered, so tiny fragments don't count as "loci"
TEMPLATE_COLS = (128, 174)  # from inputs/TR_domains.all_tr_seqs.tsv
PAD = 20


def all_hits(tbl_path):
    hits = []
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            target, hmmfrom, hmmto, alifrom, alito, strand, evalue = (
                f[0], int(f[4]), int(f[5]), int(f[6]), int(f[7]), f[11], float(f[12]),
            )
            if evalue <= EVALUE_THRESHOLD and (hmmto - hmmfrom) >= MIN_SPAN:
                hits.append(dict(target=target, hmmfrom=hmmfrom, hmmto=hmmto,
                                  alifrom=alifrom, alito=alito, strand=strand, evalue=evalue))
    return hits


def parse_stockholm_rf_and_seqs(sto_path):
    seqs = defaultdict(str)
    rf = ""
    with open(sto_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GC RF"):
                rf += line.split(None, 2)[2]
            elif line.startswith("#=GR") or line.startswith("#=GC") or line.startswith("#") or not line.strip() or line == "//":
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    seqs[parts[0]] += parts[1]
    return rf, seqs


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    species, genome_path, tbl_path, hmm_path, outdir = sys.argv[1:]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{species}.tsv"

    hits = all_hits(tbl_path)
    if len(hits) < 2:
        out_path.write_text(f"{species}\tTOO_FEW_LOCI\t{len(hits)}\t\t\t\t\n")
        return

    tmpdir = Path(tempfile.mkdtemp(prefix=f"trloci_{species}_"))
    try:
        fa = tmpdir / "genome.fa"
        if genome_path.endswith(".gz"):
            with gzip.open(genome_path, "rt") as fin, open(fa, "w") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copyfile(genome_path, fa)
        subprocess.run([ESL_SFETCH, "--index", str(fa)], check=True, capture_output=True)

        gdf_lines = []
        for i, h in enumerate(hits):
            lo, hi = sorted([h["alifrom"], h["alito"]])
            g_start, g_end = lo - PAD, hi + PAD
            if h["strand"] == "-":
                g_start, g_end = hi + PAD, lo - PAD
            gdf_lines.append(f"locus{i}\t{g_start}\t{g_end}\t{h['target']}\n")
        gdf = tmpdir / "loci.gdf"
        gdf.write_text("".join(gdf_lines))
        loci_fa = tmpdir / "loci.fa"
        subprocess.run([ESL_SFETCH, "-Cf", "-o", str(loci_fa), str(fa), str(gdf)], check=True, capture_output=True)

        sto = tmpdir / "loci.sto"
        subprocess.run([HMMALIGN, "--dna", "--trim", "-o", str(sto), hmm_path, str(loci_fa)],
                        check=True, capture_output=True)
        rf, seqs = parse_stockholm_rf_and_seqs(str(sto))
        match_positions = [i for i, c in enumerate(rf) if c == "x"]

        def template_cols(seq):
            idxs = match_positions[TEMPLATE_COLS[0] - 1:TEMPLATE_COLS[1]]
            return "".join(seq[i] for i in idxs if i < len(seq))

        rows = []
        for i, h in enumerate(hits):
            key = f"locus{i}"
            aligned = seqs.get(key, "")
            template_aligned = template_cols(aligned)
            template_ungapped = template_aligned.replace("-", "").upper()
            rows.append((species, "OK", h["target"], f"{h['alifrom']}-{h['alito']}", h["strand"],
                         f"{h['evalue']:.2e}", template_ungapped))

        with open(out_path, "w") as out:
            for r in rows:
                out.write("\t".join(r) + "\n")
        print(f"[extract_tr_template_all_loci] {species}: {len(rows)} loci -> {out_path}", file=sys.stderr)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
