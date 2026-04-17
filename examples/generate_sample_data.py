#!/usr/bin/env python3
"""
generate_sample_data.py
-----------------------
Generates a synthetic daily stock-return panel for testing csestudy.

DGP:
    r_{it} = alpha_i + beta_i' F_t + epsilon_{it} + abnormal_{it}

    F_t ~ N(mu_F, Sigma_F)          3 factors (market, size, value)
    beta_i ~ N(beta_bar, diag(s^2)) firm-specific loadings
    epsilon_{it} ~ N(0, sigma_i^2)  idiosyncratic, sigma_i ~ LogNormal
    abnormal_{it} = 0.15% * z(LNMV_i) on event day 0 (size-dependent)

Trading dates are drawn from actual S&P 500 trading days (market.stbcal),
centered around 2008-10-06 as event day 0.

Panel structure mirrors CRSP conventions:
    permno   : firm identifier (10001, 10002, ...)
    date     : S&P 500 trading dates from market.stbcal
    ret      : daily return (decimal, e.g. 0.012 = 1.2%)
    prc      : simulated closing price (for abs(prc) > 5 filter)
    sprtrn   : S&P 500 proxy return (market factor)
    lag_LNMV : lagged log market value (cross-sectional covariate)
    sic3     : 3-digit SIC code stub
    treated  : 1 = treated firm, 0 = control

Calibration targets (approximate):
    Avg firm daily vol        ~2.0%
    Cross-sectional ret corr  ~0.25
    Market factor daily vol   ~1.1%
    Idiosyncratic share       ~70% of variance

Output:
    examples/sample_data.csv
    examples/sample_data.dta

Usage:
    python generate_sample_data.py [--seed 20260416] [--nfirms 300] [--csv-only]

Example csestudy invocation (Stata):
    use examples/sample_data.dta, clear
    bcal create trading, from(date) gen(trading_date) center(20081006) replace
    tsset permno trading_date
    csestudy ret lag_LNMV if abs(prc) > 5, ///
        eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)
"""

import argparse
from datetime import date
import numpy as np
import pandas as pd
from pathlib import Path


def _patch_stata_format(dta_path: Path, varname: str,
                        columns: list[str], fmt: str) -> None:
    """Patch a variable's display format in a Stata 118 .dta file.

    pandas to_stata always writes int32 as 'long' with format '%12.0g'.
    This overwrites that format string in-place (e.g., to '%td') so Stata
    displays the column correctly without needing recast or format commands.
    """
    import struct
    col_idx = columns.index(varname)
    with open(dta_path, "r+b") as f:
        raw = f.read()
        # .dta 118: format strings are 57 bytes each, in <formats> block
        fmt_start = raw.find(b"<formats>") + len(b"<formats>")
        offset = fmt_start + col_idx * 57
        encoded = fmt.encode("ascii").ljust(57, b"\x00")
        f.seek(offset)
        f.write(encoded)


# ── S&P 500 Trading Calendar (from market.stbcal) ───────────────────────
# Range: 2006-01-03 to 2009-12-31
# Omits weekends + market holidays
_HOLIDAYS = [
    # 2006
    "2006-01-16", "2006-02-20", "2006-04-14", "2006-05-29", "2006-07-04",
    "2006-09-04", "2006-11-23", "2006-12-25",
    # 2007
    "2007-01-01", "2007-01-02", "2007-01-15", "2007-02-19", "2007-04-06",
    "2007-05-28", "2007-07-04", "2007-09-03", "2007-11-22", "2007-12-25",
    # 2008
    "2008-01-01", "2008-01-21", "2008-02-18", "2008-03-21", "2008-05-26",
    "2008-07-04", "2008-09-01", "2008-11-27", "2008-12-25",
    # 2009
    "2009-01-01", "2009-01-19", "2009-02-16", "2009-04-10", "2009-05-25",
    "2009-07-03", "2009-09-07", "2009-11-26", "2009-12-25",
]


def _sp500_trading_dates() -> list[date]:
    """Return sorted list of S&P 500 trading dates, 2006-01-03 to 2009-12-31."""
    holidays = set(pd.to_datetime(_HOLIDAYS).date)
    bdays = pd.bdate_range("2006-01-03", "2009-12-31")
    return sorted(d.date() for d in bdays if d.date() not in holidays)


def generate_sample_data(
    seed: int = 20260416,
    n_firms: int = 300,
    pre_event_days: int = 440,
    post_event_days: int = 20,
    event_date: date = date(2008, 10, 6),
    abnormal_return: float = 0.0015,
) -> pd.DataFrame:
    """Generate synthetic CRSP-like panel data.

    Parameters
    ----------
    seed : int
        RNG seed for reproducibility.
    n_firms : int
        Number of firms in the cross-section.
    pre_event_days : int
        Trading days before event day 0.
    post_event_days : int
        Trading days after event day 0.
    event_date : date
        Calendar date of event day 0.
    abnormal_return : float
        True abnormal return injected on event day (decimal).

    Returns
    -------
    pd.DataFrame
        Panel with columns: permno, date, ret, prc, sprtrn, lag_LNMV,
        sic3, treated.
    """
    rng = np.random.default_rng(seed)

    # ── Trading calendar from stbcal ─────────────────────────────────────
    all_dates = _sp500_trading_dates()
    event_idx = all_dates.index(event_date)
    start_idx = event_idx - pre_event_days
    end_idx = event_idx + post_event_days
    assert start_idx >= 0, f"Not enough pre-event dates (need {pre_event_days})"
    assert end_idx < len(all_dates), f"Not enough post-event dates"
    dates = all_dates[start_idx : end_idx + 1]
    n_days = len(dates)
    event_day_index = pre_event_days  # 0-based index of event day

    # ── Factor structure ─────────────────────────────────────────────────
    # 3 factors: market, size, value
    # Daily means (annualized: ~8%, ~3%, ~3%)
    mu_f = np.array([0.0003, 0.00012, 0.00012])

    # Daily covariance (market vol ~1.1%, size ~0.5%, value ~0.5%)
    daily_vol = np.array([0.011, 0.005, 0.005])
    corr_f = np.array([
        [1.00, 0.20, 0.15],
        [0.20, 1.00, 0.30],
        [0.15, 0.30, 1.00],
    ])
    cov_f = np.outer(daily_vol, daily_vol) * corr_f

    # Draw factor returns: (n_days, 3)
    factors = rng.multivariate_normal(mu_f, cov_f, size=n_days)

    # S&P proxy = market factor (column 0)
    sprtrn = factors[:, 0]

    # ── Firm characteristics ─────────────────────────────────────────────
    # Log market value: N(8.5, 1.5) in log-millions → range ~$10M to $100B
    # Drawn first so factor loadings can depend on size
    lnmv = rng.normal(8.5, 1.5, size=n_firms)

    # Factor loadings — correlated with size (realistic: small firms have
    # higher market beta, load more on size/value factors). This correlation
    # is what makes GLS outperform OLS: common factor shocks create spurious
    # day-to-day variation in the OLS coefficient on lag_LNMV that GLS removes.
    lnmv_centered = lnmv - lnmv.mean()
    beta_market = 1.4 - 0.05 * lnmv + rng.normal(0, 0.20, size=n_firms)
    beta_size = 0.6 - 0.04 * lnmv_centered + rng.normal(0, 0.20, size=n_firms)
    beta_value = 0.3 - 0.02 * lnmv_centered + rng.normal(0, 0.20, size=n_firms)
    betas = np.column_stack([beta_market, beta_size, beta_value])

    # Idiosyncratic volatility: log-normal for cross-sectional dispersion
    # Median ~1.6% daily, mean ~1.8%, some firms up to 4%
    sigma_idio = np.exp(rng.normal(np.log(0.016), 0.35, size=n_firms))

    # Small positive daily alpha (most near zero)
    alpha = rng.normal(0.0001, 0.0003, size=n_firms)

    # Initial price: drawn from realistic distribution, most $10–$200
    # Log-normal centered around $40
    prc_initial = np.exp(rng.normal(np.log(40), 0.7, size=n_firms))

    # SIC codes: non-financial industries
    sic3_options = np.array([
        201, 283, 357, 366, 367, 382, 384, 481, 489,
        737, 738, 739, 131, 211, 262, 281, 308, 331,
        341, 355, 371, 421, 451, 531, 541, 581, 701,
        781, 871, 874,
    ])
    sic3 = rng.choice(sic3_options, size=n_firms)

    # Firm identifiers
    permnos = np.arange(10001, 10001 + n_firms)

    # Treatment indicator: above-median size firms get positive AR,
    # below-median get negative AR (the event is size-dependent)
    treated = (lnmv >= np.median(lnmv)).astype(int)

    # ── Generate returns ─────────────────────────────────────────────────
    systematic = betas @ factors.T          # (n_firms, n_days)
    idio = rng.normal(0, 1, size=(n_firms, n_days)) * sigma_idio[:, None]
    ret = alpha[:, None] + systematic + idio

    # Inject abnormal return on event day: size-dependent for all firms.
    # The abnormal return varies with lag_LNMV so csestudy detects an
    # unusual cross-sectional gradient on the event day (the typical use
    # case: an event with differential effects by firm characteristic).
    # AR_i = abnormal_return * (lnmv_i - mean(lnmv)) / std(lnmv)
    lnmv_z = (lnmv - lnmv.mean()) / lnmv.std()
    ret[:, event_day_index] += abnormal_return * lnmv_z

    # ── Simulate price paths (for abs(prc) > 5 filter) ──────────────────
    # Cumulate returns from initial price; prc can go negative (CRSP
    # convention: negative prc = bid-ask midpoint, but abs() is used)
    prc = np.zeros((n_firms, n_days))
    prc[:, 0] = prc_initial * (1 + ret[:, 0])
    for t in range(1, n_days):
        prc[:, t] = prc[:, t - 1] * (1 + ret[:, t])

    # A handful of penny stocks: set ~5% of firms to low initial price
    n_penny = int(0.05 * n_firms)
    penny_idx = rng.choice(n_firms, size=n_penny, replace=False)
    for i in penny_idx:
        prc[i, :] = prc[i, :] / prc[i, :].mean() * rng.uniform(1.5, 4.5)

    # ── Assemble panel ───────────────────────────────────────────────────
    # Stata epoch: 01jan1960
    _STATA_EPOCH = date(1960, 1, 1)
    stata_dates = np.array([(d - _STATA_EPOCH).days for d in dates],
                           dtype=np.int32)

    rows = []
    for i in range(n_firms):
        firm_df = pd.DataFrame({
            "permno": np.int32(permnos[i]),
            "date": stata_dates,
            "ret": np.round(ret[i, :], 6),
            "prc": np.round(prc[i, :], 2),
            "sprtrn": np.round(sprtrn, 6),
            "lag_LNMV": round(lnmv[i], 4),
            "sic3": np.int16(sic3[i]),
            "treated": np.int8(treated[i]),
        })
        rows.append(firm_df)

    panel = pd.concat(rows, ignore_index=True)

    return panel


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic CRSP-like panel for csestudy testing."
    )
    parser.add_argument("--seed", type=int, default=20260416)
    parser.add_argument("--nfirms", type=int, default=300)
    parser.add_argument("--csv-only", action="store_true",
                        help="Skip .dta output (avoids pyreadstat dependency)")
    args = parser.parse_args()

    print(f"Generating panel: {args.nfirms} firms, seed={args.seed}")
    panel = generate_sample_data(seed=args.seed, n_firms=args.nfirms)

    outdir = Path(__file__).parent
    csv_path = outdir / "sample_data.csv"
    panel.to_csv(csv_path, index=False)
    print(f"Saved {csv_path}  ({len(panel):,} obs)")

    if not args.csv_only:
        try:
            dta_panel = panel.copy()
            # date is already Stata int32 (days since 01jan1960).
            # Write without convert_dates so it stays as long, then
            # patch the %td format directly in the .dta binary.
            dta_path = outdir / "sample_data.dta"
            dta_panel.to_stata(
                dta_path,
                write_index=False,
                version=118,  # Stata 14+
                variable_labels={
                    "permno": "Firm identifier",
                    "date": "Trading date",
                    "ret": "Daily return (decimal)",
                    "prc": "Closing price (synthetic)",
                    "sprtrn": "Market index return (decimal)",
                    "lag_LNMV": "Lagged log market value",
                    "sic3": "3-digit SIC code",
                    "treated": "Event treatment indicator",
                },
            )
            # Patch the date column's display format from %12.0g to %td.
            # The .dta 118 format stores formats as 57-byte fixed-width
            # strings in a block starting after the typlist and varnames.
            _patch_stata_format(dta_path, "date",
                                dta_panel.columns.tolist(), "%td")
            print(f"Saved {dta_path}")
        except Exception as e:
            print(f"Skipping .dta (pandas.to_stata failed: {e})")
            print("Run with --csv-only or install pyreadstat.")

    # ── Summary statistics ───────────────────────────────────────────────
    event_date = date(2008, 10, 6)
    print(f"\n── Summary ────────────────────────────────────────────")
    print(f"  Firms:          {panel['permno'].nunique()}")
    print(f"  Trading days:   {panel['date'].nunique()}")
    print(f"  Date range:     {panel['date'].min()} to {panel['date'].max()}")
    print(f"  Observations:   {len(panel):,}")
    print(f"  Above-med size: {panel.loc[panel['treated']==1, 'permno'].nunique()}")
    print(f"  Event date:     {event_date}")
    print(f"  Mean daily ret: {panel['ret'].mean()*100:.3f}%")
    print(f"  Std daily ret:  {panel['ret'].std()*100:.3f}%")
    print(f"  Penny stocks:   {(panel.groupby('permno')['prc'].apply(lambda x: x.abs().mean() < 5)).sum()}")
    cross_ret = panel.pivot(index="date", columns="permno", values="ret")
    avg_corr = cross_ret.corr().values
    np.fill_diagonal(avg_corr, np.nan)
    print(f"  Avg pairwise ρ: {np.nanmean(avg_corr):.3f}")
    print(f"  Injected AR:    ±0.15% × z(lag_LNMV) on event day (size-dependent)")
    print(f"\n  Stata usage:")
    print(f"    use examples/sample_data.dta, clear")
    print(f"    bcal create trading, from(date) gen(trading_date) center(20081006) replace")
    print(f"    tsset permno trading_date")
    print(f"    csestudy ret lag_LNMV if abs(prc)>5, ///")
    print(f"        eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)")


if __name__ == "__main__":
    main()
