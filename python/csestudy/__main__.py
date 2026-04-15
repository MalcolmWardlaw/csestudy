"""
CLI entry point for csestudy.

Usage:
    python -m csestudy --csv data.csv --event-date 0 --pre-start -200 --pre-end -1 \
        --depvar ret --indepvars lag_LNMV --panel-var permno --time-var date \
        --method gls --npc 100
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd

from csestudy.core import CSEventStudy


def _try_int(val: str):
    """Attempt to parse a string as an integer; return string otherwise."""
    try:
        return int(val)
    except ValueError:
        return val


def main():
    p = argparse.ArgumentParser(
        prog="csestudy",
        description="Time-series inference for cross-sectional event studies",
    )
    p.add_argument("--csv", required=True, help="Path to CSV input file")
    p.add_argument("--event-date", required=True, help="Event date (integer or date string)")
    p.add_argument("--pre-start", required=True, help="First pre-event date")
    p.add_argument("--pre-end", required=True, help="Last pre-event date")
    p.add_argument("--depvar", default="ret", help="Dependent variable column (default: ret)")
    p.add_argument("--indepvars", nargs="*", default=[], help="Independent variable columns")
    p.add_argument("--panel-var", default="permno", help="Panel identifier column (default: permno)")
    p.add_argument("--time-var", default="date", help="Time variable column (default: date)")
    p.add_argument("--sample", default=None, help="Pandas query string for subsetting")
    p.add_argument("--method", choices=["ols", "gls"], default="ols", help="Estimation method")
    p.add_argument("--npc", type=int, default=100, help="Number of principal components (GLS)")
    p.add_argument("--solver", choices=["cholesky", "woodbury"], default="cholesky",
                    help="GLS solver (default: cholesky)")
    p.add_argument("--out-prefix", default=None, help="Output file prefix (writes JSON + CSV)")
    args = p.parse_args()

    print(f"[info] Loading {args.csv}", file=sys.stderr)
    df = pd.read_csv(args.csv)

    model = CSEventStudy(
        df,
        event_date=_try_int(args.event_date),
        pre_start=_try_int(args.pre_start),
        pre_end=_try_int(args.pre_end),
        depvar=args.depvar,
        indepvars=args.indepvars if args.indepvars else None,
        panel_var=args.panel_var,
        time_var=args.time_var,
        sample=args.sample,
        method=args.method,
        npc=args.npc,
        solver=args.solver,
    )

    result = model.fit(verbose=True)
    print(result.summary())

    if args.out_prefix:
        summary = {
            "method": result.method,
            "params": {n: float(b) for n, b in zip(result.param_names, result.params)},
            "p_cdf": {n: float(p) for n, p in zip(result.param_names, result.p_cdf)},
            "p_parametric": {n: float(p) for n, p in zip(result.param_names, result.p_parametric)},
            "n_obs_event": result.n_obs_event,
            "n_pre_event_days": result.n_pre_event_days,
        }
        json_path = f"{args.out_prefix}_summary.json"
        csv_path = f"{args.out_prefix}_pre_betas.csv"
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)
        beta_df = pd.DataFrame(result.pre_event_betas, columns=result.param_names)
        beta_df.to_csv(csv_path, index=False)
        print(f"[done] Summary: {json_path}", file=sys.stderr)
        print(f"[done] Pre-betas: {csv_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
