#!/usr/bin/env python3
"""
Overlay the per-sample P(s) contact-decay curves on one log-log figure.
A reference slope of -1 is drawn for orientation (healthy mammalian Hi-C).
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaling", nargs="+", required=True, help="per-sample scaling TSVs")
    ap.add_argument("--names", required=True, help="comma-separated sample names (same order)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = args.names.split(",")
    if len(names) != len(args.scaling):
        raise SystemExit("[ERROR] --names count must match --scaling count")

    fig, ax = plt.subplots(figsize=(6, 5))
    ycol = None
    for name, path in zip(names, args.scaling):
        df = pd.read_csv(path, sep="\t")
        ycol = "balanced_avg" if "balanced_avg" in df.columns else "count_avg"
        d = df[(df["s_bp"] > 0) & (df[ycol] > 0)]
        ax.plot(d["s_bp"], d[ycol], label=name, lw=1.2)

    # reference P(s) ~ s^-1 anchored to the middle of the x-range
    if ycol is not None:
        import numpy as np
        xs = np.array(ax.get_xlim())
        x0 = xs[0] if xs[0] > 0 else 1e4
        ax.plot([x0, x0 * 100], [1, 0.01], "k--", lw=0.8, alpha=0.6, label="slope -1")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("genomic distance s (bp)")
    ax.set_ylabel("contact frequency P(s) (balanced, normalised)")
    ax.set_title("Hi-C contact-distance decay")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"[OK] wrote {args.out}")


if __name__ == "__main__":
    main()
