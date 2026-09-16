#!/usr/bin/env python3
"""
Compartmentalisation STRENGTH from a saddle plot (cooltools API).

A saddle stratifies observed/expected cis contacts by the compartment
eigenvalue (E1) of the two interacting bins. Strong compartmentalisation =>
high A-A and B-B (corners) and low A-B (off-diagonal). We report:

    strength = mean(A-A, B-B corners) / mean(A-B, B-A corners)

using the most-extreme `--extreme-fraction` of quantile bins for the corners.
Also saves the saddle matrix (.npz) for plotting.
"""
import argparse
import numpy as np
import pandas as pd
import cooler
import cooltools
from cooltools.api.saddle import saddle


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcool", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--eigs", required=True, help="cooltools eigs-cis vecs.tsv (must contain E1)")
    ap.add_argument("--n-bins", type=int, default=38)
    ap.add_argument("--extreme-fraction", type=float, default=0.2,
                    help="fraction of quantile bins forming each saddle corner")
    ap.add_argument("--view", default="",
                    help="optional BED4 viewframe (chrom start end name); "
                         "default = whole chromosomes carrying E1")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--out-strength", required=True)
    args = ap.parse_args()

    clr = cooler.Cooler(f"{args.mcool}::resolutions/{args.resolution}")

    eig = pd.read_csv(args.eigs, sep="\t")
    eig = eig[["chrom", "start", "end", "E1"]]

    # viewframe: from --view if given, else whole chromosomes that carry E1
    if args.view:
        view_df = pd.read_csv(args.view, sep="\t", header=None,
                              names=["chrom", "start", "end", "name"])
        view_df = view_df[view_df["chrom"].isin(set(eig["chrom"]))].reset_index(drop=True)
    else:
        chroms = [c for c in clr.chromnames if c in set(eig["chrom"])]
        view_df = pd.DataFrame({
            "chrom": chroms,
            "start": 0,
            "end": [clr.chromsizes[c] for c in chroms],
            "name": chroms,
        })

    expected = cooltools.expected_cis(clr, view_df=view_df, nproc=1)

    q = (0.025, 0.975)
    interaction_sum, interaction_count = saddle(
        clr, expected, eig, "cis",
        n_bins=args.n_bins, qrange=q, view_df=view_df,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        sad = interaction_sum / interaction_count

    np.savez(args.out_prefix + ".saddledump.npz",
             saddle=sad, sum=interaction_sum, count=interaction_count)

    # corners: ignore the first/last row+col (outlier bins added by cooltools)
    core = sad[1:-1, 1:-1]
    n = core.shape[0]
    k = max(1, int(round(n * args.extreme_fraction)))
    AA = np.nanmean(core[:k, :k])        # low E1 x low E1  (B-B)
    BB = np.nanmean(core[-k:, -k:])      # high E1 x high E1 (A-A)
    AB = np.nanmean(core[:k, -k:])
    BA = np.nanmean(core[-k:, :k])
    strength = (AA + BB) / (AB + BA)

    with open(args.out_strength, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"resolution\t{args.resolution}\n")
        fh.write(f"extreme_fraction\t{args.extreme_fraction}\n")
        fh.write(f"AA_BB_mean\t{(AA + BB) / 2:.6g}\n")
        fh.write(f"AB_BA_mean\t{(AB + BA) / 2:.6g}\n")
        fh.write(f"compartment_strength\t{strength:.6g}\n")
    print(f"[OK] compartment strength = {strength:.4g}")


if __name__ == "__main__":
    main()
