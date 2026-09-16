#!/usr/bin/env python3
"""loop_caller_agreement.py

Compute slop-tolerant 2-D agreement between chromatin-loop callers
(chromosight, Mustache, HiCCUPS) for one sample.

A loop is an anchor pair (chrA, posA) <-> (chrB, posB). Two callers "agree"
on a loop when the pairs (chrA, posA, chrB, posB) match within +/- slop_bp
on BOTH anchors AND on the same chromosome pair.

Outputs:
    --out_tsv     :  one-row-per-sample TSV with per-caller counts +
                    pairwise Jaccards + consensus counts at >=2-of-N
                    and >=3-of-N support.
    --out_bedpe   :  BEDPE5 of consensus loops (chr1,start1,end1,
                    chr2,start2,end2,callers,support).
"""

from __future__ import annotations
import argparse, csv, gzip, os, sys
from collections import defaultdict
from itertools import combinations


# ---------------------------------------------------------------------------
# Parsers — each caller has a slightly different output schema
# ---------------------------------------------------------------------------
def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_chromosight(path):
    """chromosight TSV: chrom1, start1, end1, chrom2, start2, end2, bin1, bin2,
    score, pvalue, qvalue, ... (header present)."""
    if not path or not os.path.isfile(path):
        return []
    pairs = []
    with _open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return []
        for row in reader:
            if len(row) < 6:
                continue
            try:
                c1, s1, e1, c2, s2, e2 = row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5])
            except (ValueError, IndexError):
                continue
            pairs.append((c1, (s1 + e1) // 2, c2, (s2 + e2) // 2))
    return pairs


def parse_mustache(path):
    """Mustache TSV: chrom1, start1, end1, chrom2, start2, end2, fdr,
    detection_scale (sometimes header, sometimes not)."""
    if not path or not os.path.isfile(path):
        return []
    pairs = []
    with _open(path) as fh:
        for raw in fh:
            row = raw.rstrip("\n").split("\t")
            if len(row) < 6:
                continue
            # header detection: first column non-numeric pos won't parse as int below
            try:
                c1, s1, e1, c2, s2, e2 = row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5])
            except (ValueError, IndexError):
                continue
            pairs.append((c1, (s1 + e1) // 2, c2, (s2 + e2) // 2))
    return pairs


def parse_hiccups(path):
    """HiCCUPS merged_loops.bedpe: chr1, x1, x2, chr2, y1, y2, ... (comment header)."""
    if not path or not os.path.isfile(path):
        return []
    pairs = []
    with _open(path) as fh:
        for raw in fh:
            if raw.startswith("#") or raw.startswith("chr1"):
                continue
            row = raw.rstrip("\n").split("\t")
            if len(row) < 6:
                continue
            try:
                c1, s1, e1, c2, s2, e2 = row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5])
            except (ValueError, IndexError):
                continue
            pairs.append((c1, (s1 + e1) // 2, c2, (s2 + e2) // 2))
    return pairs


# ---------------------------------------------------------------------------
# Slop-tolerant matching between two pair-sets
# ---------------------------------------------------------------------------
def _bin(pos, slop):
    # Place into a bucket of size 2*slop so two points within slop fall in
    # adjacent buckets (we'll check the neighbour bucket too).
    return pos // max(1, 2 * slop)


def jaccard_2d(A, B, slop):
    """Jaccard of two loop sets under +/- slop on BOTH anchors."""
    if not A or not B:
        return 0.0
    # Build an index of B by (chrom1, chrom2, bin1, bin2)
    idx = defaultdict(list)
    for (c1, p1, c2, p2) in B:
        idx[(c1, c2, _bin(p1, slop), _bin(p2, slop))].append((p1, p2))
    matched_in_B = set()
    matches_A = 0
    for (c1, pA1, c2, pA2) in A:
        b1, b2 = _bin(pA1, slop), _bin(pA2, slop)
        hit = False
        for db1 in (-1, 0, 1):
            for db2 in (-1, 0, 1):
                for (pB1, pB2) in idx.get((c1, c2, b1 + db1, b2 + db2), ()):
                    if abs(pA1 - pB1) <= slop and abs(pA2 - pB2) <= slop:
                        hit = True
                        matched_in_B.add((c1, pB1, c2, pB2))
        if hit:
            matches_A += 1
    inter = (matches_A + len(matched_in_B)) / 2.0
    union = len(A) + len(B) - inter
    return inter / union if union > 0 else 0.0


# ---------------------------------------------------------------------------
# Multi-caller consensus via single-linkage clustering
# ---------------------------------------------------------------------------
def consensus(callers, slop):
    """
    callers: dict {name -> [(c1, p1, c2, p2), ...]}
    Returns list of (chr1, pos1, chr2, pos2, [caller_names], support_count)
    """
    # Pool all loops with their source caller
    pool = []
    for name, pairs in callers.items():
        for (c1, p1, c2, p2) in pairs:
            pool.append((c1, p1, c2, p2, name))

    # Sort by chromosome pair + bins so contiguous neighbours appear close
    pool.sort(key=lambda x: (x[0], x[2], x[1], x[3]))

    # Single-linkage cluster by sliding through and merging anything within slop
    parent = list(range(len(pool)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    # Index pool by chrom-pair + 2D-bin for nearest-neighbour merging
    idx = defaultdict(list)
    for k, (c1, p1, c2, p2, name) in enumerate(pool):
        idx[(c1, c2, _bin(p1, slop), _bin(p2, slop))].append(k)
    for (c1, c2, b1, b2), members in idx.items():
        for db1 in (-1, 0, 1):
            for db2 in (-1, 0, 1):
                neigh = idx.get((c1, c2, b1 + db1, b2 + db2), ())
                for ki in members:
                    p1i, p2i = pool[ki][1], pool[ki][3]
                    for kj in neigh:
                        if kj <= ki:
                            continue
                        p1j, p2j = pool[kj][1], pool[kj][3]
                        if abs(p1i - p1j) <= slop and abs(p2i - p2j) <= slop:
                            union(ki, kj)

    # Collect cluster members
    clusters = defaultdict(list)
    for k in range(len(pool)):
        clusters[find(k)].append(pool[k])

    result = []
    for members in clusters.values():
        caller_names = sorted(set(m[4] for m in members))
        c1 = members[0][0]
        c2 = members[0][2]
        p1 = sum(m[1] for m in members) // len(members)
        p2 = sum(m[3] for m in members) // len(members)
        result.append((c1, p1, c2, p2, caller_names, len(caller_names)))
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--slop_bp", type=int, default=10000)
    p.add_argument("--chromosight", default="")
    p.add_argument("--mustache",    default="")
    p.add_argument("--hiccups",     default="")
    p.add_argument("--out_tsv",     required=True)
    p.add_argument("--out_bedpe",   required=True)
    args = p.parse_args()

    pairs = {}
    if args.chromosight and os.path.isfile(args.chromosight):
        pairs["chromosight"] = parse_chromosight(args.chromosight)
    if args.mustache and os.path.isfile(args.mustache):
        pairs["mustache"] = parse_mustache(args.mustache)
    if args.hiccups and os.path.isfile(args.hiccups):
        pairs["hiccups"] = parse_hiccups(args.hiccups)

    if not pairs:
        sys.exit(f"[loop_caller_agreement] no caller outputs found for sample={args.sample}")

    callers = sorted(pairs.keys())
    counts = {c: len(pairs[c]) for c in callers}

    # Pairwise Jaccards
    jac = {}
    for a, b in combinations(callers, 2):
        jac[f"jaccard_{a}_{b}"] = jaccard_2d(pairs[a], pairs[b], args.slop_bp)

    # Consensus clusters
    clusters = consensus(pairs, args.slop_bp)
    consensus_counts = defaultdict(int)
    for c in clusters:
        consensus_counts[c[5]] += 1   # support count
    n_callers = len(callers)
    cons_cols = {}
    for k in range(1, n_callers + 1):
        cons_cols[f"n_consensus_{k}of{n_callers}"] = consensus_counts.get(k, 0)

    # Fraction of each caller's calls that landed in >=2-of-N consensus
    high_conf_set = set()
    for (c1, p1, c2, p2, names, support) in clusters:
        if support >= 2:
            for n in names:
                high_conf_set.add((n, c1, p1, c2, p2))
    frac_cols = {}
    for c in callers:
        in_consensus = sum(1 for (nm, _, _, _, _) in high_conf_set if nm == c)
        frac_cols[f"frac_{c}_in_consensus"] = (in_consensus / counts[c]) if counts[c] else 0.0

    # Write per-sample row TSV
    cols = ["sample", "slop_bp", "n_callers"] + \
           [f"n_{c}" for c in callers] + \
           list(jac.keys()) + \
           list(cons_cols.keys()) + \
           list(frac_cols.keys())
    vals = [args.sample, args.slop_bp, n_callers] + \
           [counts[c] for c in callers] + \
           [f"{jac[k]:.4f}" for k in jac] + \
           [cons_cols[k] for k in cons_cols] + \
           [f"{frac_cols[k]:.4f}" for k in frac_cols]

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(str(v) for v in vals) + "\n")

    # Write consensus BEDPE (>=2-of-N)
    with open(args.out_bedpe, "w") as fh:
        fh.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tcallers\tsupport\n")
        slop = args.slop_bp
        for (c1, p1, c2, p2, names, support) in clusters:
            if support < 2:
                continue
            # Represent as +/- slop intervals around each anchor midpoint.
            fh.write(f"{c1}\t{max(0, p1 - slop // 2)}\t{p1 + slop // 2}\t"
                     f"{c2}\t{max(0, p2 - slop // 2)}\t{p2 + slop // 2}\t"
                     f"{','.join(names)}\t{support}\n")

    sys.stderr.write(f"[loop_caller_agreement] sample={args.sample} "
                     f"n_callers={n_callers} consensus_2of{n_callers}={cons_cols.get(f'n_consensus_2of{n_callers}', 0)}\n")


if __name__ == "__main__":
    main()
