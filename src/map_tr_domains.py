#!/usr/bin/env python3
"""Map plant TR gene domains (USE, TATA-box, G-rich 5' end, Template,
Conserved region / template boundary element, C-rich 3' end, Terminator -
per Fajkus et al. 2019 NAR, https://doi.org/10.1093/nar/gkz695, Fig 1B/1D)
onto TR.hmm match-state columns (1..LENG), so they can be sliced out of any
genome hit scored against that HMM using its reported hmmfrom/hmmto.

Method:
  1. Template is anchored directly: for a set of reference species with a
     known telomere repeat unit (Table S2 / Fig 1B of the paper), find the
     longest exact match in their TR sequence to any rotation/register of
     revcomp(repeat), length in [len(repeat)+1, 2*len(repeat)].
  2. hmmalign the 75-sequence reference set (inputs/TR.fasta) against
     inputs/TR.hmm, giving every reference sequence's positions in the
     HMM's 1..LENG match-column coordinate system (Stockholm convention:
     match column = uppercase letter or '-'; insert column = lowercase or
     '.'). This locates the Template anchors in HMM-column space.
  3. Build a per-column conservation profile (mode-base frequency, GC
     content) across all 75 references, smooth it, and call contiguous
     high-conservation blocks. Blocks are labelled by their position
     relative to the anchored Template and by composition (AT-rich near
     the 5' end -> TATA-box; GC-rich just before Template -> G-rich 5'
     end; conserved block right after Template -> Conserved region /
     template boundary element; C-rich block further 3' -> C-rich 3' end;
     near-pure T-run at the modeled 3' edge -> Terminator; the single
     near-invariant block at the very 5' edge -> USE).

Output: inputs/TR_domains.tsv (domain, hmm_col_start, hmm_col_end,
consensus, mean_conservation, mean_gc, evidence).

This is a heuristic structural map, not a base-pair-exact annotation -
Template is empirically anchored and high-confidence; Terminator (poly-T
run) is essentially unambiguous; USE/TATA-box/G-rich/Conserved
region/C-rich are conservation+composition calls and should be treated as
approximate windows, wide enough to contain the true element.
"""
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_ROOT / "inputs"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
# The reference 75-species domain layout (Fajkus et al. 2019) is defined
# against inputs/TR.hmm (LENG 317), the original seed HMM. But every
# per-species genome scan in outputs/tr_tbls/*_all_tr_seqs_020226.tbl was
# run against outputs/tr_tbls/all_tr_seqs.hmm (LENG 346) - the HMM rebuilt
# from ~859 candidate hits after the first identification round (see
# run_tr_hmms.bash). Domain columns must be mapped onto whichever HMM the
# hit coordinates you're using were actually scored against - pass it as
# argv[1], defaulting to the production all_tr_seqs.hmm since that's what
# extract_tr_domains.py consumes.
TR_HMM = Path(sys.argv[1]) if len(sys.argv) > 1 else OUTPUT_DIR / "tr_tbls" / "all_tr_seqs.hmm"
TR_FASTA = INPUT_DIR / "TR.fasta"
HMMALIGN = "/software/team301/hmmer-3.4/src/hmmalign"
OUT_TSV = Path(sys.argv[2]) if len(sys.argv) > 2 else INPUT_DIR / f"TR_domains.{TR_HMM.stem}.tsv"

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# species substring (must appear in a TR.fasta header) -> known telomere
# repeat unit, from Fajkus et al. 2019 Table S2 / Fig 1B and standard plant
# telomere literature (canonical TTTAGGG for Arabidopsis/Nicotiana/Posidonia).
ANCHOR_REPEATS = {
    "Allium angulosum": "CTCGGTTATGGG",
    "Allium cepa": "CTCGGTTATGGG",
    "Allium ericetorum": "CTCGGTTATGGG",
    "Allium fistulosum": "CTCGGTTATGGG",
    "Allium nutans": "CTCGGTTATGGG",
    "Allium ursinum": "CTCGGTTATGGG",
    "Cestrum elegans": "TTTTTTAGGG",
    "Scilla peruviana": "TTAGGG",
    "Tulbaghia violacea": "TTAGGG",
    "Arabidopsis thaliana": "TTTAGGG",
    "Arabidopsis lyrata": "TTTAGGG",
    "Nicotiana sylvestris": "TTTAGGG",
    "Posidonia oceanica": "TTTAGGG",
}

def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def rotations(s):
    return {s[i:] + s[:i] for i in range(len(s))}


def load_fasta(path):
    recs = {}
    hdr, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr:
                    recs[hdr] = "".join(seq)
                hdr, seq = line[1:], []
            else:
                seq.append(line)
        if hdr:
            recs[hdr] = "".join(seq)
    return recs


def best_repeat_match(seq, repeat):
    rc = revcomp(repeat)
    best = None
    for rot in rotations(rc):
        doubled = rot + rot
        for length in range(len(repeat) + 1, 2 * len(repeat) + 1):
            cand = doubled[:length]
            for m in re.finditer(re.escape(cand), seq):
                if best is None or length > best[0]:
                    best = (length, m.start(), m.end())
    return best


def find_template_anchors():
    recs = load_fasta(TR_FASTA)
    anchors = {}
    for sp, repeat in ANCHOR_REPEATS.items():
        hdr = next((h for h in recs if sp in h), None)
        if hdr is None:
            print(f"[map_tr_domains] WARNING: anchor species {sp!r} not found in TR.fasta", file=sys.stderr)
            continue
        hit = best_repeat_match(recs[hdr], repeat)
        if hit is None:
            print(f"[map_tr_domains] WARNING: no template match for {sp}", file=sys.stderr)
            continue
        _, start, end = hit
        anchors[sp] = (hdr, start, end)
    return anchors


def run_hmmalign():
    with tempfile.NamedTemporaryFile(suffix=".sto", delete=False) as tmp:
        sto_path = tmp.name
    subprocess.run(
        [HMMALIGN, "--dna", "--trim", "-o", sto_path, str(TR_HMM), str(TR_FASTA)],
        check=True,
    )
    return sto_path


def parse_stockholm(sto_path):
    id_to_desc = {}
    seqs = defaultdict(str)
    with open(sto_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GS") and " DE " in line:
                m = re.match(r"#=GS\s+(\S+)\s+DE\s+(.*)", line)
                if m:
                    id_to_desc[m.group(1)] = m.group(2)
            elif line and not line.startswith("#") and line != "//":
                parts = line.split(None, 1)
                if len(parts) == 2:
                    sid, chunk = parts
                    seqs[sid] += chunk
    return id_to_desc, seqs


def hmm_length(hmm_path):
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith("LENG"):
                return int(line.split()[1])
    raise ValueError("LENG not found in HMM file")


def ungapped_to_matchcol(aln, ustart, uend):
    """Map an ungapped [ustart, uend) residue range to HMM match columns."""
    match_col = 0
    ungapped_pos = 0
    start_col = end_col = None
    for ch in aln:
        is_match_col = ch.isupper() or ch == "-"
        if is_match_col:
            match_col += 1
        if ch not in "-.":
            if ustart <= ungapped_pos < uend and is_match_col:
                if start_col is None:
                    start_col = match_col
                end_col = match_col
            ungapped_pos += 1
    return start_col, end_col


def build_match_matrix(seqs, hmm_len):
    matrix = []
    for sid, aln in seqs.items():
        col_bases = [None] * hmm_len
        match_col = 0
        for ch in aln:
            is_match_col = ch.isupper() or ch == "-"
            if is_match_col:
                match_col += 1
                if ch.isupper():
                    col_bases[match_col - 1] = ch
        matrix.append(col_bases)
    return matrix


def conservation_profile(matrix, hmm_len):
    profile = []
    for col in range(hmm_len):
        bases = [row[col] for row in matrix if row[col] is not None]
        n = len(bases)
        if n == 0:
            profile.append((0, "-", 0.0, 0.0))
            continue
        counts = Counter(bases)
        top_base, top_count = counts.most_common(1)[0]
        frac = top_count / n
        gc = (counts.get("G", 0) + counts.get("C", 0)) / n
        profile.append((n, top_base, frac, gc))
    return profile


def smooth(values, window=3):
    n = len(values)
    out = []
    half = window // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        out.append(sum(values[lo:hi]) / (hi - lo))
    return out


def call_blocks(frac, threshold=0.65, min_len=4, max_gap=2):
    """Contiguous runs of smoothed conservation >= threshold, 0-based cols,
    merging runs separated by small gaps."""
    n = len(frac)
    raw = [i for i in range(n) if frac[i] >= threshold]
    if not raw:
        return []
    blocks = []
    start = prev = raw[0]
    for i in raw[1:]:
        if i - prev <= max_gap:
            prev = i
        else:
            blocks.append((start, prev))
            start = prev = i
    blocks.append((start, prev))
    return [(s, e) for s, e in blocks if e - s + 1 >= min_len]


def main():
    print("[map_tr_domains] finding template anchors by known repeat unit...", file=sys.stderr)
    anchors = find_template_anchors()
    for sp, (hdr, s, e) in anchors.items():
        print(f"  {sp:25s} template@{s}-{e}", file=sys.stderr)

    print("[map_tr_domains] running hmmalign...", file=sys.stderr)
    sto_path = run_hmmalign()
    id_to_desc, seqs = parse_stockholm(sto_path)
    desc_to_id = {d: s for s, d in id_to_desc.items()}
    hmm_len = hmm_length(TR_HMM)

    raw_cols = {}
    for sp, (hdr, ustart, uend) in anchors.items():
        sid = next((s for d, s in desc_to_id.items() if sp in d), None)
        if sid is None:
            continue
        start_col, end_col = ungapped_to_matchcol(seqs[sid], ustart, uend)
        print(f"  {sp:25s} HMM cols {start_col}-{end_col}", file=sys.stderr)
        if start_col is not None:
            raw_cols[sp] = (start_col, end_col)

    # Cluster anchors by their column midpoint and keep only the largest
    # tight cluster. This HMM may be built from real genome hits (see
    # below) that only cover some reference lineages well - e.g. a
    # production HMM trained on DToL species will align exotic
    # divergent-repeat references (Allium, Cestrum, Scilla...) to a
    # different, mutually-consistent-but-wrong part of the model, forming
    # a second cluster rather than scattering as simple outliers. Taking
    # the largest contiguous cluster (consecutive-midpoint gap <=
    # max_gap) is robust to that bimodality without hardcoding species.
    max_gap = 20
    midpoints = sorted((((s + e) / 2), sp) for sp, (s, e) in raw_cols.items())
    clusters, current = [], [midpoints[0]]
    for prev, cur in zip(midpoints, midpoints[1:]):
        if cur[0] - prev[0] <= max_gap:
            current.append(cur)
        else:
            clusters.append(current)
            current = [cur]
    clusters.append(current)
    best_cluster = max(clusters, key=len)
    kept = {sp: raw_cols[sp] for _, sp in best_cluster}
    dropped = set(raw_cols) - set(kept)
    if dropped:
        print(f"[map_tr_domains] dropped (separate, less-supported cluster; {len(kept)} anchors kept): "
              f"{', '.join(sorted(dropped))}", file=sys.stderr)
    template_cols = list(kept.values())

    tmpl_start = min(c[0] for c in template_cols)
    tmpl_end = max(c[1] for c in template_cols)
    print(f"[map_tr_domains] consensus Template window (robust anchors): {tmpl_start}-{tmpl_end}", file=sys.stderr)

    matrix = build_match_matrix(seqs, hmm_len)
    profile = conservation_profile(matrix, hmm_len)
    frac = [p[2] for p in profile]
    gc = [p[3] for p in profile]
    consensus_base = [p[1] for p in profile]
    frac_smooth = smooth(frac, window=3)

    blocks = call_blocks(frac_smooth, threshold=0.65, min_len=4, max_gap=2)
    print(f"[map_tr_domains] {len(blocks)} conserved blocks called (smoothed frac >= 0.65):", file=sys.stderr)
    for s, e in blocks:
        cons = "".join(consensus_base[s:e + 1])
        mean_gc = sum(gc[s:e + 1]) / (e - s + 1)
        mean_frac = sum(frac[s:e + 1]) / (e - s + 1)
        print(f"    cols {s+1}-{e+1}  cons={cons}  mean_frac={mean_frac:.2f}  mean_gc={mean_gc:.2f}", file=sys.stderr)

    # Classify blocks by position relative to Template + composition.
    domains = {}
    domains["Template"] = (tmpl_start, tmpl_end,
                            f"anchored: revcomp(known telomere repeat) exact match, "
                            f"{len(template_cols)}/{len(raw_cols)} reference species after outlier rejection")

    pre_template = [(s, e) for s, e in blocks if e + 1 < tmpl_start]
    post_template = [(s, e) for s, e in blocks if s + 1 > tmpl_end]

    if pre_template:
        use_block = pre_template[0]
        domains["USE"] = (use_block[0] + 1, use_block[1] + 1,
                           "conservation: near-invariant block at 5'-most edge of alignment")
        tata_candidates = [b for b in pre_template[1:] if sum(1 for c in "AT" if consensus_base[b[0]] in "AT")]
    tata_candidates = [b for b in pre_template[1:]] if len(pre_template) > 1 else []
    # TATA-box: AT-rich block among the pre-template conserved blocks (excluding USE)
    def at_frac(block):
        s, e = block
        return sum(1 for i in range(s, e + 1) if consensus_base[i] in "AT") / (e - s + 1)

    tata_candidates = sorted(tata_candidates, key=at_frac, reverse=True)
    if tata_candidates and at_frac(tata_candidates[0]) >= 0.6:
        b = tata_candidates[0]
        domains["TATA-box"] = (b[0] + 1, b[1] + 1, "composition: AT-rich conserved block, downstream of USE")

    # G-rich 5' end: GC-rich pre-template block closest to Template
    grich_candidates = sorted(
        [b for b in pre_template if b not in tata_candidates[:1] and b != pre_template[0]],
        key=lambda b: -b[1],
    )
    grich_candidates = [b for b in grich_candidates if sum(gc[b[0]:b[1] + 1]) / (b[1] - b[0] + 1) >= 0.6]
    if grich_candidates:
        b = grich_candidates[0]
        domains["G-rich 5' end"] = (b[0] + 1, b[1] + 1, "composition: GC-rich conserved block immediately preceding Template")

    if post_template:
        cr_block = post_template[0]
        domains["Conserved region"] = (cr_block[0] + 1, cr_block[1] + 1,
                                        "conservation: conserved block immediately following Template (candidate template boundary element)")

    # Terminator: block(s) at the very 3' edge of the alignment, near-pure T
    terminator_candidates = [b for b in blocks if b[1] >= hmm_len - 12]
    terminator_candidates = [b for b in terminator_candidates
                              if sum(1 for i in range(b[0], b[1] + 1) if consensus_base[i] == "T") / (b[1] - b[0] + 1) >= 0.8]
    if terminator_candidates:
        b = terminator_candidates[-1]
        domains["Terminator"] = (b[0] + 1, b[1] + 1, "composition: near-pure poly-T run at the 3'-most edge of the alignment (classic Pol III terminator)")

    # C-rich 3' end: C-rich block between Conserved region and Terminator
    term_start = domains.get("Terminator", (hmm_len + 1,))[0]
    cr_end = domains.get("Conserved region", (0, 0))[1]
    crich_candidates = [b for b in post_template if b[0] + 1 > cr_end and b[1] + 1 < term_start]
    def c_frac(block):
        s, e = block
        return sum(1 for i in range(s, e + 1) if consensus_base[i] == "C") / (e - s + 1)
    crich_candidates = sorted(crich_candidates, key=c_frac, reverse=True)
    if crich_candidates and c_frac(crich_candidates[0]) >= 0.5:
        b = crich_candidates[0]
        domains["C-rich 3' end"] = (b[0] + 1, b[1] + 1, "composition: C-rich conserved block between Conserved region and Terminator")

    order = ["USE", "TATA-box", "G-rich 5' end", "Template", "Conserved region", "C-rich 3' end", "Terminator"]
    with open(OUT_TSV, "w") as out:
        out.write("domain\thmm_col_start\thmm_col_end\tconsensus\tevidence\n")
        for name in order:
            if name not in domains:
                print(f"[map_tr_domains] WARNING: could not confidently call domain {name!r}", file=sys.stderr)
                continue
            s, e, evidence = domains[name]
            cons = "".join(consensus_base[s - 1:e])
            out.write(f"{name}\t{s}\t{e}\t{cons}\t{evidence}\n")

    print(f"[map_tr_domains] wrote {OUT_TSV}", file=sys.stderr)
    Path(sto_path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
