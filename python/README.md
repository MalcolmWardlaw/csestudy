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

df = pd.read_csv("my_panel.csv")

# OLS
result = CSEventStudy(
    df,
    event_date=0,
    pre_start=-200,
    pre_end=-1,
    depvar="ret",
    indepvars=["lag_LNMV"],
    panel_var="permno",
    time_var="date",
    sample="abs_prc > 5",      # optional pandas query (like Stata's [if])
    method="ols",
).fit()

print(result.summary())
print(result.params)       # coefficient vector
print(result.p_cdf)        # empirical CDF p-values

# GLS with Cholesky (default, most precise)
result_gls = CSEventStudy(
    df, event_date=0, pre_start=-200, pre_end=-1,
    depvar="ret", indepvars=["lag_LNMV"],
    method="gls", npc=100,
).fit()

# GLS with Woodbury identity (faster, slightly less precise)
result_wb = CSEventStudy(
    df, event_date=0, pre_start=-200, pre_end=-1,
    depvar="ret", indepvars=["lag_LNMV"],
    method="gls", npc=100, solver="woodbury",
).fit()
```

## Command Line

```bash
python -m csestudy \
    --csv my_panel.csv \
    --event-date 0 --pre-start -200 --pre-end -1 \
    --depvar ret --indepvars lag_LNMV \
    --method gls --npc 100 \
    --out-prefix results
```

## Dependencies

- Python >= 3.11
- numpy >= 1.24
- pandas >= 2.0
- scipy >= 1.11
