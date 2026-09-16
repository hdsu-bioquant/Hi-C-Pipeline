#!/usr/bin/env python3
"""
Collate the key library metrics from pairtools dedup `.stats.txt` files into one
cross-sample table for a quick comparison and for inclusion in reports.

Columns:
    sample  total  mapped  valid_pairs  dup_rate  cis_frac  long_cis_frac(>=20kb)  trans_frac
"""
import argparse
import pandas as pd


def parse_stats(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, val = line.split("\t", 1)
            d[key] = val
    return d


def num(d, key, default=0.0):
    try:
        return float(d.get(key, default))
    except ValueError:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", nargs="+", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = args.names.split(",")
    if len(names) != len(args.stats):
        raise SystemExit("[ERROR] --names count must match --stats count")

    rows = []
    for name, path in zip(names, args.stats):
        d = parse_stats(path)
        total   = num(d, "total")
        mapped  = num(d, "total_mapped")
        dups    = num(d, "total_dups")
        nodups  = num(d, "total_nodups")
        cis     = num(d, "cis")
        trans   = num(d, "trans")
        # pairtools writes cumulative long-cis bins like cis_20kb+
        long_cis = num(d, "cis_20kb+")
        cis_trans = cis + trans
        rows.append({
            "sample": name,
            "total": int(total),
            "mapped": int(mapped),
            "valid_pairs": int(nodups),
            "dup_rate": round(dups / mapped, 4) if mapped else float("nan"),
            "cis_frac": round(cis / cis_trans, 4) if cis_trans else float("nan"),
            "long_cis_frac_ge20kb": round(long_cis / cis, 4) if cis else float("nan"),
            "trans_frac": round(trans / cis_trans, 4) if cis_trans else float("nan"),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"[OK] wrote {args.out}")


if __name__ == "__main__":
    main()
