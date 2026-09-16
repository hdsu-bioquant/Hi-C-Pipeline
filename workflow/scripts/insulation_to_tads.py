#!/usr/bin/env python3
"""
Convert a cooltools `insulation` table into TAD interval and boundary BED files.

cooltools writes, per genomic bin and per diamond window W, the columns:
    log2_insulation_score_<W>, boundary_strength_<W>, is_boundary_<W>
A bin flagged `is_boundary_<W> == True` is a TAD boundary (a local insulation
minimum that passed the strength threshold). TADs are the segments between
consecutive boundaries (chromosome ends act as virtual boundaries).
"""
import argparse
import re
import pandas as pd


def pick(df, prefix, win):
    exact = f"{prefix}_{win}"
    if exact in df.columns:
        return exact
    cand = [c for c in df.columns if c.startswith(prefix + "_") and str(win) in c]
    if cand:
        return cand[0]
    raise KeyError(f"No column like {prefix}_{win} in {list(df.columns)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--insulation", required=True)
    ap.add_argument("--window-bp", type=int, required=True)
    ap.add_argument("--keep-chroms", default="",
                    help="comma-separated chromosomes to keep (default: all)")
    ap.add_argument("--out-tads", required=True)
    ap.add_argument("--out-boundaries", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.insulation, sep="\t", comment="#")
    need = {"chrom", "start", "end"}
    if not need.issubset(df.columns):
        raise SystemExit(f"[ERROR] insulation table missing {need}; has {list(df.columns)}")

    is_bound_col = pick(df, "is_boundary", args.window_bp)
    try:
        strength_col = pick(df, "boundary_strength", args.window_bp)
    except KeyError:
        strength_col = None

    # drop bad bins (cooltools marks low-coverage/blacklisted bins)
    bad = df["is_bad_bin"].astype(bool) if "is_bad_bin" in df.columns else pd.Series(False, index=df.index)
    df = df.loc[~bad].copy()

    if args.keep_chroms:
        keep = set(args.keep_chroms.split(","))
        df = df[df["chrom"].isin(keep)].copy()
    if df.empty:
        raise SystemExit("[ERROR] no bins left after filtering")

    # ---- boundaries BED ----
    is_b = df[is_bound_col].astype(bool)
    bdf = df.loc[is_b, ["chrom", "start", "end"]].copy()
    if strength_col:
        bdf["strength"] = df.loc[is_b, strength_col].values
    bdf.to_csv(args.out_boundaries, sep="\t", header=False, index=False)

    # ---- TAD intervals (segments between consecutive boundary starts) ----
    rows = []
    for chrom, g in df.groupby("chrom", sort=False):
        g = g.sort_values("start")
        cuts = [int(g["start"].iloc[0])]
        cuts += [int(x) for x in g.loc[g[is_bound_col].astype(bool), "start"].tolist()]
        cuts.append(int(g["end"].iloc[-1]))
        cuts = sorted(set(cuts))
        for a, b in zip(cuts[:-1], cuts[1:]):
            if b > a:
                rows.append((chrom, a, b))
    pd.DataFrame(rows, columns=["chrom", "start", "end"]).to_csv(
        args.out_tads, sep="\t", header=False, index=False
    )
    print(f"[OK] {len(rows)} TADs, {len(bdf)} boundaries (window={args.window_bp})")


if __name__ == "__main__":
    main()
