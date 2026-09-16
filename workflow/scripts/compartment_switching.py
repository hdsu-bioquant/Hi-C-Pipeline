#!/usr/bin/env python3
"""
Per-bin A/B compartment switching across the time-course.

Merges the per-sample E1 eigenvector tracks on (chrom,start,end), records the
compartment label (A if E1>0 else B) for each sample, and flags bins that
change label at least once. Output is a tidy table you can filter for
switching bins or feed into downstream enrichment.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vecs", nargs="+", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = args.names.split(",")
    if len(names) != len(args.vecs):
        raise SystemExit("[ERROR] --names must match --vecs count")

    merged = None
    for name, path in zip(names, args.vecs):
        df = pd.read_csv(path, sep="\t")[["chrom", "start", "end", "E1"]]
        df = df.rename(columns={"E1": f"E1_{name}"})
        merged = df if merged is None else merged.merge(df, on=["chrom", "start", "end"], how="outer")

    e1_cols = [f"E1_{n}" for n in names]
    lab_cols = []
    for n in names:
        lab = f"comp_{n}"
        merged[lab] = np.where(merged[f"E1_{n}"] > 0, "A",
                       np.where(merged[f"E1_{n}"] < 0, "B", np.nan))
        lab_cols.append(lab)

    def n_states(row):
        vals = {row[c] for c in lab_cols if isinstance(row[c], str)}
        return len(vals)

    merged["n_compartment_states"] = merged.apply(n_states, axis=1)
    merged["switches"] = (merged["n_compartment_states"] > 1).astype(int)

    merged = merged.sort_values(["chrom", "start"])
    merged.to_csv(args.out, sep="\t", index=False)
    n_sw = int(merged["switches"].sum())
    n_tot = int((merged["n_compartment_states"] > 0).sum())
    print(f"[OK] {n_sw}/{n_tot} assignable bins switch A<->B across {len(names)} timepoints")


if __name__ == "__main__":
    main()
