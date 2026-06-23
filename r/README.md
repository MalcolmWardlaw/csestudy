# csestudy (R)

R implementation of the time-series approach to cross-sectional event study inference.

> Cohn, Johnson, Liu, and Wardlaw (2026), "Past is Prologue: Inference from the Cross Section of Returns Around an Event," *Journal of Financial Economics* 180, 104278.

This package mirrors the Stata (`csestudy`) and Python (`CSEventStudy`) implementations in the same repository and reproduces their numerical output on the shared sample data.

## Installation

```r
# from GitHub (parallel to  pip install ./python)
remotes::install_github("MalcolmWardlaw/csestudy", subdir = "r")

# or from a local clone
# R CMD INSTALL r
```

The package has no hard dependencies beyond base R and `stats`.

## Quick Start

```r
library(csestudy)

df <- read.csv("examples/sample_data.csv")

# Reindex calendar dates to sequential trading-day integers centered on the
# event date (the R analogue of Stata's `bcal ... center(...)`). The event,
# pre_start, and pre_end values below refer to this sequential index.
tds <- sort(unique(df$date))
df$tdate <- match(df$date, tds) - match(17811L, tds)   # 17811 = 2008-10-06 -> 0

# OLS with time-series-corrected errors
fit <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                event_date = 0, pre_start = -200, pre_end = -1,
                panel_var = "permno", time_var = "tdate",
                sample = abs(prc) > 5,        # optional row filter, like Stata's [if]
                method = "ols")

summary(fit)        # coefficient / CDF p-value / parametric p-value table
coef(fit)           # event-date coefficient vector
fit$p_cdf           # empirical CDF p-values

# GLS with 100 principal components (Cholesky, default, most precise)
fit_gls <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                    event_date = 0, pre_start = -200, pre_end = -1,
                    panel_var = "permno", time_var = "tdate",
                    sample = abs(prc) > 5, method = "gls", npc = 100)

# GLS with the Woodbury identity (faster, slightly less precise)
fit_wb <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                   event_date = 0, pre_start = -200, pre_end = -1,
                   panel_var = "permno", time_var = "tdate",
                   sample = abs(prc) > 5, method = "gls", npc = 100,
                   solver = "woodbury")
```

## Arguments

| Argument | Description |
|----------|-------------|
| `data` | Long-form panel data frame (one row per panel × time). |
| `depvar`, `indepvars` | Dependent-variable column and a character vector of regressors (`NULL` for intercept only). |
| `event_date`, `pre_start`, `pre_end` | Values of `time_var` for the event date and the pre-event window. Require `event_date > pre_end > pre_start`. |
| `panel_var`, `time_var` | Panel-id and time columns. `time_var` should be a sequential integer with non-trading days omitted. |
| `sample` | Optional unquoted expression to subset rows before estimation (e.g. `abs(prc) > 5`). |
| `method` | `"ols"` (default) or `"gls"`. |
| `npc` | Number of principal components for the GLS covariance matrix (default 100). |
| `solver` | `"cholesky"` (default) or `"woodbury"`. Requires `method = "gls"`. |

See `?csestudy` for full details and the package help for the GLS data requirements (a strongly balanced panel of non-missing returns across each pre-event window).

## Dependencies

- R >= 3.5
- `stats` (base R)
