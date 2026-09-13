#!/usr/bin/env python3
"""For one species, extract each TR gene domain (USE, TATA-box, G-rich 5'
end, Template, Conserved region, C-rich 3' end, Terminator - see
map_tr_domains.py / inputs/TR_domains.tsv) from its genome assembly, using
the best nhmmer TR.hmm hit's hmmfrom/hmmto -> alifrom/alito mapping to
interpolate genome coordinates for each domain's HMM match-column range.

For the Template domain specifically, also cross-check against that
species' empirically observed telomeric repeat (outputs/tidk/*.tidk.tsv,
already discovered directly from the genome by TIDK, independent of TR
gene homology): does the padded, extracted window contain the reverse
complement of (a rotation of) the TIDK repeat?

Usage:
  extract_tr_domains.py <species> <genome_fa_gz> <nhmmer_tbl> <domains_tsv> [tidk_tsv] <outdir>

Writes <outdir>/<species>.tsv (one row per domain found).
"""
import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ESL_SFETCH = "/software/team301/hmmer-3.4/easel/miniapps/esl-sfetch"
PAD = 15  # bp of padding on each side of the interpolated domain window
EVALUE_THRESHOLD = 1e-4

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def rotations(s):
    return {s[i:] + s[:i] for i in range(len(s))}


def load_domains(path):
    domains = []
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            name, start, end, cons, evidence = line.rstrip("\n").split("\t")
            domains.append((name, int(start), int(end)))
    return domains


def best_hit(tbl_path):
    """Lowest-E-value row from an nhmmer --tblout file, or None."""
    best = None
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            target, hmmfrom, hmmto, alifrom, alito, strand, evalue = (
                f[0], int(f[4]), int(f[5]), int(f[6]), int(f[7]), f[11], float(f[12]),
            )
            if evalue > EVALUE_THRESHOLD:
                continue
            if best is None or evalue < best["evalue"]:
                best = dict(target=target, hmmfrom=hmmfrom, hmmto=hmmto,
                            alifrom=alifrom, alito=alito, strand=strand, evalue=evalue)
    return best


def interpolate(hit, hmm_col):
    """Genome coordinate (1-based) for a given HMM match-column, by linear
    interpolation across the hit's aligned span."""
    span_hmm = hit["hmmto"] - hit["hmmfrom"]
    if span_hmm == 0:
        return hit["alifrom"]
    frac = (hmm_col - hit["hmmfrom"]) / span_hmm
    return round(hit["alifrom"] + frac * (hit["alito"] - hit["alifrom"]))


TIDK_TOP_N = 15  # TIDK ranks repeats genome-wide by total count, so the
# most abundant one is often a non-telomeric satellite rather than the true
# (comparatively low-copy, chromosome-end-restricted) telomere repeat -
# check the top N rows, not just rank 1.


def load_tidk_repeats(path, top_n=TIDK_TOP_N):
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        next(fh)  # header
        return [line.split("\t")[0].strip() for line in fh.readlines()[:top_n] if line.strip()]


def find_repeat_in_window(seq, repeat):
    """Does seq contain the reverse complement of any rotation of repeat?
    Returns the matched substring and its span, or None."""
    rc = revcomp(repeat)
    best = None
    for rot in rotations(rc):
        doubled = rot + rot
        for length in range(len(repeat), 2 * len(repeat) + 1):
            cand = doubled[:length]
            idx = seq.find(cand)
            if idx >= 0 and (best is None or length > best[0]):
                best = (length, idx, cand)
    return best


def main():
    if len(sys.argv) not in (6, 7):
        sys.exit(__doc__)
    if len(sys.argv) == 6:
        species, genome_path, tbl_path, domains_path, outdir = sys.argv[1:]
        tidk_path = None
    else:
        species, genome_path, tbl_path, domains_path, tidk_path, outdir = sys.argv[1:]

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{species}.tsv"

    domains = load_domains(domains_path)
    hit = best_hit(tbl_path)

    rows = []
    if hit is None:
        out_path.write_text(f"{species}\tNO_CONFIDENT_HIT\t\t\t\t\t\t\n")
        print(f"[extract_tr_domains] {species}: no hit passing E<={EVALUE_THRESHOLD}", file=sys.stderr)
        return

    # Check the canonical land-plant repeat first (most probable hypothesis,
    # regardless of what TIDK ranked highest), then TIDK's candidates in
    # its own rank order (rank-1 by genome-wide count is not always the
    # telomere repeat - satellite DNA can outnumber it, hence checking
    # more than just rank 1). Take the FIRST candidate (in this priority
    # order) with any match, not the longest match across all candidates -
    # a short, low-information repeat (e.g. 4-5bp) can "coincidentally"
    # produce a longer doubled-rotation match than a real 7bp canonical
    # hit purely because it has more short rotations to try, which would
    # wrongly out-rank a genuine canonical match.
    tidk_repeats = load_tidk_repeats(tidk_path)
    candidate_repeats = list(dict.fromkeys(["TTTAGGG"] + tidk_repeats))

    tmpdir = Path(tempfile.mkdtemp(prefix=f"trdom_{species}_"))
    try:
        fa = tmpdir / "genome.fa"
        if genome_path.endswith(".gz"):
            with gzip.open(genome_path, "rt") as fin, open(fa, "w") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            # esl-sfetch --index writes a .ssi next to the fasta; the
            # source data directory is typically read-only, so copy to
            # tmpdir first rather than indexing in place.
            shutil.copyfile(genome_path, fa)
        subprocess.run([ESL_SFETCH, "--index", str(fa)], check=True, capture_output=True)

        gdf_lines = []
        windows = {}
        for name, col_start, col_end in domains:
            if hit["hmmfrom"] > col_start or hit["hmmto"] < col_end:
                continue  # hit doesn't span this domain's columns
            g_start = interpolate(hit, col_start) - PAD
            g_end = interpolate(hit, col_end) + PAD
            if hit["strand"] == "-":
                g_start, g_end = interpolate(hit, col_start) + PAD, interpolate(hit, col_end) - PAD
            key = re.sub(r"[^A-Za-z0-9]+", "_", name)
            windows[key] = name
            gdf_lines.append(f"{key}\t{g_start}\t{g_end}\t{hit['target']}\n")

        if not gdf_lines:
            out_path.write_text(f"{species}\tHIT_TOO_PARTIAL\t{hit['target']}\t{hit['hmmfrom']}-{hit['hmmto']}\t\t\t\t\n")
            return

        gdf = tmpdir / "windows.gdf"
        gdf.write_text("".join(gdf_lines))
        extracted = subprocess.run(
            [ESL_SFETCH, "-Cf", str(fa), str(gdf)], check=True, capture_output=True, text=True
        ).stdout

        seqs = {}
        cur_name, cur_seq = None, []
        for line in extracted.splitlines():
            if line.startswith(">"):
                if cur_name:
                    seqs[cur_name] = "".join(cur_seq)
                cur_name, cur_seq = line[1:].split()[0], []
            else:
                cur_seq.append(line.strip())
        if cur_name:
            seqs[cur_name] = "".join(cur_seq)

        for key, name in windows.items():
            seq = seqs.get(key, "")
            tidk_note = ""
            if name == "Template" and candidate_repeats:
                match, repeat_hit, src = None, None, None
                for repeat in candidate_repeats:
                    m = find_repeat_in_window(seq, repeat)
                    if m:
                        match, repeat_hit = m, repeat
                        if repeat == "TTTAGGG":
                            src = "canonical TTTAGGG"
                        else:
                            src = f"tidk rank {tidk_repeats.index(repeat) + 1}"
                        break
                if match:
                    tidk_note = f"revcomp({repeat_hit}) MATCH len={match[0]} [{src}]"
                else:
                    tidk_note = f"no match in top {len(tidk_repeats)} tidk repeats or canonical TTTAGGG"
            rows.append((species, "OK", hit["target"], f"{hit['hmmfrom']}-{hit['hmmto']}",
                         name, seq, str(len(seq)), tidk_note))

        with open(out_path, "w") as out:
            for r in rows:
                out.write("\t".join(r) + "\n")
        print(f"[extract_tr_domains] {species}: wrote {len(rows)} domain(s) to {out_path}", file=sys.stderr)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
