#!/usr/bin/env python3
"""
Build a boundary-presence matrix across timepoints.

Boundaries from different samples rarely sit on the exact same bp, so we cluster
boundaries that fall within `--merge-window` of each other into a consensus
position, then record which samples contribute to each consensus boundary.
Output columns: chrom  start  end  <sample1> ... <sampleN>  n_samples
(1 = boundary present in that sample, 0 = absent).
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beds", nargs="+", required=True, help="boundary BEDs, same order as --names")
    ap.add_argument("--names", required=True, help="comma-separated sample names")
    ap.add_argument("--merge-window", type=int, default=10000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = args.names.split(",")
    if len(names) != len(args.beds):
        raise SystemExit("[ERROR] --names count must match --beds count")

    frames = []
    for name, bed in zip(names, args.beds):
        df = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
        df["sample"] = name
        frames.append(df)
    allb = pd.concat(frames, ignore_index=True).sort_values(["chrom", "start"])

    rows = []
    for chrom, g in allb.groupby("chrom", sort=False):
        g = g.sort_values("start")
        cluster = []
        anchor = None
        for r in g.itertuples(index=False):
            if anchor is None or r.start - anchor <= args.merge_window:
                cluster.append(r)
                anchor = r.start if anchor is None else anchor
            else:
                rows.append(_emit(chrom, cluster, names))
                cluster = [r]
                anchor = r.start
        if cluster:
            rows.append(_emit(chrom, cluster, names))

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] {len(out)} consensus boundaries across {len(names)} samples")


def _emit(chrom, cluster, names):
    starts = [c.start for c in cluster]
    present = {c.sample for c in cluster}
    rec = {"chrom": chrom, "start": min(starts), "end": max(c.end for c in cluster)}
    for n in names:
        rec[n] = int(n in present)
    rec["n_samples"] = sum(rec[n] for n in names)
    return rec


if __name__ == "__main__":
    main()
