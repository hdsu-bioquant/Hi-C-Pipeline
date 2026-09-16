#!/usr/bin/env python3
"""
Per-sample boundary agreement across TAD callers.

Reads boundary sets from up to three callers (cooltools insulation, Juicer
Arrowhead, HiCExplorer hicFindTADs) for ONE sample, computes pairwise Jaccard-
like agreement under a slop tolerance, and emits:

  - <sample>.agreement.tsv : one row, one column per metric (n boundaries per
                             caller, pairwise Jaccards, consensus-support
                             counts, fraction of each caller's boundaries in
                             the >=2-of-N consensus).
  - <sample>.consensus.bed : merged boundary intervals called by >= 2 callers.
                             BED5: chrom, start, end, callers (comma list),
                             score = number of supporting callers (2 or 3).

Boundary point extraction:
  - From a BOUNDARY BED (insulation, hicFindTADs) -> each row's midpoint.
  - From a DOMAIN BED (Arrowhead)                  -> each row's start AND end.

Agreement metric:
  Single-linkage clustering across callers on the same chromosome within
  ±slop_bp. For each cluster, the supporting-caller count is len(unique
  callers in the cluster). Jaccard(A,B) is computed symmetrically as
  mean(|A overlapped by B|, |B overlapped by A|) / (|A| + |B| - that_mean).

Adaptive mode:
  - hicexplorer flag absent -> two-caller mode (insulation, arrowhead). The
    TSV columns and consensus counts adapt automatically (e.g.
    n_consensus_2of2 / n_consensus_1of2).
"""
import argparse
import bisect
from collections import defaultdict


# ----- BED parsing --------------------------------------------------------- #

def _read_intervals(path, keep_chroms=None):
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue  # header row or malformed
            chrom = parts[0]
            if keep_chroms is not None and chrom not in keep_chroms:
                continue
            out.append((chrom, s, e))
    return out


def boundary_points_from_boundary_bed(path, keep_chroms=None):
    """Each row is a boundary INTERVAL -> point = midpoint."""
    return [(c, (s + e) // 2) for c, s, e in _read_intervals(path, keep_chroms)]


def boundary_points_from_domain_bed(path, keep_chroms=None):
    """Each row is a DOMAIN -> two boundary points: start and end."""
    pts = []
    for c, s, e in _read_intervals(path, keep_chroms):
        pts.append((c, s))
        pts.append((c, e))
    return pts


# ----- agreement primitives ------------------------------------------------ #

def _by_chrom(points):
    bc = defaultdict(list)
    for c, p in points:
        bc[c].append(p)
    for c in bc:
        bc[c].sort()
    return bc


def overlap_count(set_a, set_b, slop):
    """How many points of A have at least one B-point within ±slop bp."""
    bbc = _by_chrom(set_b)
    n = 0
    for c, p in set_a:
        positions = bbc.get(c)
        if not positions:
            continue
        idx = bisect.bisect_left(positions, p - slop)
        if idx < len(positions) and positions[idx] <= p + slop:
            n += 1
    return n


def jaccard(set_a, set_b, slop):
    """Slop-tolerant Jaccard. Symmetric intersect = mean(A->B, B->A)."""
    if not set_a or not set_b:
        return 0.0
    ab = overlap_count(set_a, set_b, slop)
    ba = overlap_count(set_b, set_a, slop)
    inter = (ab + ba) / 2
    union = len(set_a) + len(set_b) - inter
    return inter / union if union > 0 else 0.0


def consensus(point_sets_by_caller, slop):
    """Single-linkage cluster ALL boundary points across callers within ±slop.
    Returns list of dicts with chrom, center, min, max, callers, n_callers."""
    tagged = []
    for caller, pts in point_sets_by_caller.items():
        for c, p in pts:
            tagged.append((c, p, caller))
    tagged.sort()

    clusters, cur = [], None
    for c, p, caller in tagged:
        if cur is None or cur['chrom'] != c or (p - cur['max']) > slop:
            if cur is not None:
                clusters.append(cur)
            cur = {'chrom': c, 'min': p, 'max': p,
                   'positions': [p], 'callers': {caller}}
        else:
            cur['max'] = p
            cur['positions'].append(p)
            cur['callers'].add(caller)
    if cur is not None:
        clusters.append(cur)
    for cl in clusters:
        cl['center'] = sum(cl['positions']) // len(cl['positions'])
        cl['n_callers'] = len(cl['callers'])
    return clusters


def write_consensus_bed(clusters, path, min_support=2):
    with open(path, 'w') as f:
        for cl in sorted(clusters, key=lambda x: (x['chrom'], x['center'])):
            if cl['n_callers'] < min_support:
                continue
            s = cl['min']
            e = max(cl['max'] + 1, s + 1)
            name = ",".join(sorted(cl['callers']))
            f.write(f"{cl['chrom']}\t{s}\t{e}\t{name}\t{cl['n_callers']}\n")


# ----- main ---------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample', required=True)
    ap.add_argument('--insulation', required=True,
                    help="cooltools boundaries BED (intervals -> midpoint)")
    ap.add_argument('--arrowhead', required=True,
                    help="Juicer Arrowhead domains BED (each row's start AND end are boundary points)")
    ap.add_argument('--hicexplorer', default=None,
                    help="HiCExplorer hicFindTADs boundaries BED (optional)")
    ap.add_argument('--slop-bp', type=int, required=True,
                    help="Tolerance (bp) for one caller's boundary to count as agreeing with another")
    ap.add_argument('--keep-chroms', default=None,
                    help="Comma-separated whitelist of chromosomes to consider (skip scaffolds). "
                         "Empty -> keep all.")
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-consensus', required=True)
    args = ap.parse_args()

    chrom_filter = set(args.keep_chroms.split(',')) if args.keep_chroms else None

    sets = {
        'insulation': boundary_points_from_boundary_bed(args.insulation, chrom_filter),
        'arrowhead':  boundary_points_from_domain_bed(args.arrowhead, chrom_filter),
    }
    if args.hicexplorer:
        sets['hicexplorer'] = boundary_points_from_boundary_bed(args.hicexplorer, chrom_filter)

    slop = args.slop_bp
    keys = list(sets.keys())
    n_callers_total = len(keys)

    # Pairwise Jaccards (preserve key order)
    pair_jacc = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            pair_jacc.append((a, b, jaccard(sets[a], sets[b], slop)))

    # Consensus clusters
    clusters = consensus(sets, slop)
    n_by_support = defaultdict(int)
    for cl in clusters:
        n_by_support[cl['n_callers']] += 1
    consensus_pts = [(cl['chrom'], cl['center'])
                     for cl in clusters if cl['n_callers'] >= 2]

    # Fraction of each caller's boundaries that landed in a >=2-of-N cluster
    frac_in_cons = {k: (overlap_count(pts, consensus_pts, slop) / len(pts) if pts else 0.0)
                    for k, pts in sets.items()}

    # ----- write per-sample TSV (one row) -----
    cols = ['sample', 'slop_bp', 'n_callers']
    cols += [f'n_{k}' for k in keys]
    cols += [f'jaccard_{a}_{b}' for (a, b, _) in pair_jacc]
    cols += [f'n_consensus_{n}of{n_callers_total}' for n in range(n_callers_total, 0, -1)]
    cols += [f'frac_{k}_in_consensus' for k in keys]

    row = {
        'sample': args.sample,
        'slop_bp': slop,
        'n_callers': n_callers_total,
    }
    for k in keys:
        row[f'n_{k}'] = len(sets[k])
    for a, b, v in pair_jacc:
        row[f'jaccard_{a}_{b}'] = f"{v:.4f}"
    for n in range(n_callers_total, 0, -1):
        row[f'n_consensus_{n}of{n_callers_total}'] = n_by_support.get(n, 0)
    for k in keys:
        row[f'frac_{k}_in_consensus'] = f"{frac_in_cons[k]:.4f}"

    with open(args.out_tsv, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        f.write('\t'.join(str(row[c]) for c in cols) + '\n')

    write_consensus_bed(clusters, args.out_consensus, min_support=2)

    # ----- stdout summary (captured in the rule log) -----
    print(f"[tad-agreement] sample={args.sample} "
          f"callers={','.join(keys)} slop={slop}bp "
          f"chrom_filter={'on' if chrom_filter else 'off'}")
    for k in keys:
        print(f"  n_{k}={len(sets[k])}")
    for a, b, v in pair_jacc:
        print(f"  jaccard({a},{b})={v:.4f}")
    for n in range(n_callers_total, 0, -1):
        print(f"  consensus_{n}of{n_callers_total}={n_by_support.get(n, 0)}")
    for k in keys:
        print(f"  frac_{k}_in_consensus={frac_in_cons[k]:.4f}")


if __name__ == '__main__':
    main()
