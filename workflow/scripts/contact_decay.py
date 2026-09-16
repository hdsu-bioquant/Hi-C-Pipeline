#!/usr/bin/env python3
"""
Contact-distance decay curve P(s) for one sample, from a balanced .mcool.

Uses cooltools.expected_cis (the diagonal-wise average of the ICE-balanced
matrix, restricted to the analysis view) and aggregates across regions into a
single genome-wide P(s) table:

    dist_bins  s_bp  n_valid  count_avg  balanced_avg

`balanced_avg` is the quantity to plot on log-log axes; its slope over
~0.1-1 Mb is ~ -1 for a healthy mammalian Hi-C library.
"""
import argparse
import pandas as pd
import cooler
import cooltools


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcool", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--view", default="")
    ap.add_argument("--ignore-diags", type=int, default=2)
    ap.add_argument("--nproc", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    clr = cooler.Cooler(f"{args.mcool}::resolutions/{args.resolution}")

    view_df = None
    if args.view:
        view_df = pd.read_csv(args.view, sep="\t", header=None,
                              names=["chrom", "start", "end", "name"])
        # keep only chromosomes present in the cooler
        view_df = view_df[view_df["chrom"].isin(set(clr.chromnames))].reset_index(drop=True)

    exp = cooltools.expected_cis(
        clr, view_df=view_df, nproc=args.nproc,
        clr_weight_name="weight", ignore_diags=args.ignore_diags,
    )

    # Column names differ slightly across cooltools versions; resolve robustly.
    def col(*cands):
        for c in cands:
            if c in exp.columns:
                return c
        return None

    c_dist = col("dist", "diag")
    c_nval = col("n_valid")
    c_csum = col("count.sum", "count_sum")
    c_bsum = col("balanced.sum", "balanced_sum")

    agg = {c_nval: "sum", c_csum: "sum"}
    if c_bsum:
        agg[c_bsum] = "sum"
    g = exp.groupby(c_dist, as_index=False).agg(agg)

    out = pd.DataFrame()
    out["dist_bins"] = g[c_dist]
    out["s_bp"] = g[c_dist] * args.resolution
    out["n_valid"] = g[c_nval]
    out["count_avg"] = g[c_csum] / g[c_nval].where(g[c_nval] > 0)
    if c_bsum:
        out["balanced_avg"] = g[c_bsum] / g[c_nval].where(g[c_nval] > 0)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote P(s) with {len(out)} distance bins -> {args.out}")


if __name__ == "__main__":
    main()
