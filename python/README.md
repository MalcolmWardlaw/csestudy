# csestudy (Python)

Python implementation of the time-series approach to cross-sectional event study inference.

> Cohn, Johnson, Liu, and Wardlaw (2026), "Past is Prologue: Inference from the Cross Section of Returns Around an Event," *Journal of Financial Economics* 180, 104278.

## Installation

```bash
pip install ./python          # from the repo root
# or
pip install -e ./python       # editable / development install
```

## Quick Start

```python
import pandas as pd
from csestudy import CSEventStudy

# Load the included sample data (or substitute your own panel)
df = pd.read_csv("examples/sample_data.csv")

# ── Reindex dates to sequential integers (Python equivalent of bcal) ──
# The time_var must contain values that event_date, pre_start, and
# pre_end refer to.  With CRSP-style data the simplest approach is to
# map calendar dates to sequential trading-day integers centered on the
# event date (analogous to Stata's  bcal create ... center(...)).
trading_dates = sorted(df["date"].unique())
event_raw = 17811                          # 2008-10-06 as Stata integer
seq_map = {d: i - trading_dates.index(event_raw) for i, d in enumerate(trading_dates)}
df["tdate"] = df["date"].map(seq_map)      # tdate: ..., -2, -1, 0, 1, 2, ...

# OLS
result = CSEventStudy(
    df,
    event_date=0,
    pre_start=-200,
    pre_end=-1,
    depvar="ret",
    indepvars=["lag_LNMV"],
    panel_var="permno",
    time_var="tdate",
    sample="prc.abs() > 5",   # optional pandas query (like Stata's [if])
    method="ols",
).fit()

print(result.summary())
print(result.params)       # coefficient vector
print(result.p_cdf)        # empirical CDF p-values

# GLS with Cholesky (default, most precise)
result_gls = CSEventStudy(
    df, event_date=0, pre_start=-200, pre_end=-1,
    depvar="ret", indepvars=["lag_LNMV"],
    panel_var="permno", time_var="tdate",
    sample="prc.abs() > 5",
    method="gls", npc=100,
).fit()

# GLS with Woodbury identity (faster, slightly less precise)
result_wb = CSEventStudy(
    df, event_date=0, pre_start=-200, pre_end=-1,
    depvar="ret", indepvars=["lag_LNMV"],
    panel_var="permno", time_var="tdate",
    sample="prc.abs() > 5",
    method="gls", npc=100, solver="woodbury",
).fit()
```

## Command Line

```bash
# NOTE: The CLI currently expects event-date / pre-start / pre-end to be
# literal values in the time-var column.  With the shipped sample data
# (Stata integer dates), pass the raw date values:
python -m csestudy \
    --csv examples/sample_data.csv \
    --event-date 17811 --pre-start 17519 --pre-end 17808 \
    --depvar ret --indepvars lag_LNMV \
    --panel-var permno --time-var date \
    --method gls --npc 100 \
    --out-prefix results
```

## Dependencies

- Python >= 3.11
- numpy >= 1.24
- pandas >= 2.0
- scipy >= 1.11
