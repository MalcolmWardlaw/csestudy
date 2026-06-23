# CSESTUDY: Efficient Inference for Cross-Sectional Event Studies

**Stata, Python, and R implementations** of the time-series approach to cross-sectional event study inference described in:

> Cohn, Johnson, Liu, and Wardlaw (2026), "Past is Prologue: Inference from the Cross Section of Returns Around an Event," *Journal of Financial Economics* 180, 104278. [doi:10.1016/j.jfineco.2026.104278](https://doi.org/10.1016/j.jfineco.2026.104278)

SSRN: https://ssrn.com/abstract=4296657

`csestudy` provides correct standard errors and p-values for cross-sectional event studies — abnormal-return regressions estimated on the cross section of stock returns around a corporate event date. It corrects for the strong cross-correlation structure in returns (which industry-clustered standard errors do not address) by benchmarking the event-day relationship against the distribution of the same relationship on pre-event days, using OLS or GLS with a PCA-based covariance estimate, and reporting an empirical-CDF p-value and a parametric p-value. The same method is available as a Stata command, a Python package, and an R package.

Feedback is welcome. Please open an issue if something appears to fail or work incorrectly.

## Implementations (Stata · Python · R)

The identical methodology is implemented in three languages. Pick whichever fits your workflow — results agree numerically across all three on the shared sample data.

| Language | Install | Entry point | Docs |
|----------|---------|-------------|------|
| **Stata** | `ssc install csestudy, replace` | `csestudy` command | [help file](csestudy.sthlp) · below |
| **Python** | `pip install ./python` | `csestudy.CSEventStudy` | [python/README.md](python/README.md) |
| **R** | `remotes::install_github("MalcolmWardlaw/csestudy", subdir = "r")` | `csestudy()` | [r/README.md](r/README.md) |

```stata
// Stata
csestudy depvar [indepvars] [if] , eventstartdate(string) firstpreeventdate(string) lastpreeventdate(string)
```

```python
# Python
from csestudy import CSEventStudy
res = CSEventStudy(df, event_date=0, pre_start=-200, pre_end=-1,
                   depvar="ret", indepvars=["lag_LNMV"], method="gls").fit()
```

```r
# R
library(csestudy)
fit <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                event_date = 0, pre_start = -200, pre_end = -1, method = "gls")
```

The remainder of this README documents the **Stata** command. See [python/README.md](python/README.md) and [r/README.md](r/README.md) for the Python and R packages.

## Basic Description

Standard event study methodologies usually fail to account for the strong cross-correlation structure in stock returns across firm characteristics. Clustering standard errors by industry does not address this problem. This package implements a time-series approach that benchmarks the event-period relationship against a distribution of the same relationship estimated on pre-event days, using either OLS or GLS (with PCA-based covariance estimation). Rejection criteria are computed as a parametric z-score and a p-value from the empirical CDF of pre-event coefficients.

## Syntax and Usage (Stata)

The data must be **tsset** by panel id and time. The time variable should be a sequential integer with non-trading days omitted. The simplest way to achieve this is via **bcal create**. See the Stata help file for details.

```stata
csestudy depvar [indepvars] [if] , eventstartdate() firstpreeventdate() lastpreeventdate()
```

### Options

| Option | Description |
|--------|-------------|
| `gls` | Use GLS estimation with PCA-based covariance matrix (default is OLS) |
| `npc(integer)` | Number of principal components for GLS covariance matrix (default 100) |
| `woodbury` | Use the Woodbury matrix identity for GLS instead of Cholesky decomposition. Faster (~50-66%) but slightly less numerically precise. Requires `gls`. |
| `coefsonly` | Report only the event-period coefficients without computing significance statistics |

### Data Requirements

Data from both the event window and the pre-event window should be loaded into Stata. The `[if]` condition applies to event and pre-event date observations, but not to the dependent variable used to construct the PCA covariance matrix under `gls`.

### GLS and Balancing the Pre-Period Data

The GLS option requires a strongly balanced panel of nonmissing values for the dependent variable across each pre-event window. It also requires data extending back before `firstpreeventdate` by a window equal to `eventstartdate` - `firstpreeventdate`. For example, with `eventstartdate(0)` and `firstpreeventdate(-200)`, you need data back to approximately t = -400. The balancing routine constructs a new balanced panel for each iteration, but if there are large gaps the comparison sample may become unrepresentative.

### Multi-Day Event Windows

`csestudy` tests a single event date per invocation. To test a multi-day event window (e.g., a five-day [0,5] CAR), construct a rolling cumulative return variable in your data before calling the command:

```stata
* Test using the last day of the window as the event date
gen ret5 = ret + l1.ret + l2.ret + l3.ret + l4.ret
csestudy ret5 lag_LNMV if abs(prc)>5, eventstartdate(4) firstpreeventdate(-199) lastpreeventdate(-1) gls npc(100)

* Test using the first day of the window as the event date
gen ret5 = ret + f1.ret + f2.ret + f3.ret + f4.ret
csestudy ret5 lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-204)
        lastpreeventdate(-5)

```

The pre-event pseudo-events automatically use the same return horizon (five-day cumulative returns centered on each pre-event date), so you only need to construct the variable once for the entire time series.

**Note:** When running separate single-day tests on consecutive days (e.g., day 0 and day 1), the pre-event window length must equal `eventstartdate` - `firstpreeventdate`. This means the `firstpreeventdate` shifts forward by one day for each subsequent event date. This is by design — each pseudo-event needs its own pre-event window of the same length.

### Event Date Input

The command accepts dates as integer values or Stata expressions evaluated at runtime, e.g.:

```stata
csestudy ret lag_LNMV, eventstartdate(100) ...
csestudy ret lag_LNMV, eventstartdate(bofd("mycal",mdy(9,19,2011))) ...
```

Note the importance of tracking trading days rather than calendar days. If your time variable is a calendar date, you can use the `bcal` command to create a trading day variable. For example:

```stata
    bcal create trading, from(date) gen(trading_date) center(20081006) replace
```

Alternatively, many users have historically used a simple sequential integer for the time variable, with non-trading days omitted.

```stata
    bysort permno (date): gen time = _n
    tsset permno time
    csestudy ret lag_LNMV, eventstartdate(265) ...
```

Note that this works fine as long as the panel is strongly balanced (i.e. all stock id observations begin at the same date and exist on all trading days). This program works fine with either approach, but users should take care to guarantee that the time variable is correctly specified and that the event date is correctly aligned with the time variable. The `bcal` approach is more robust to missing data and non-trading days, while the sequential integer approach can be more convenient if your data is already structured that way.



## Installation (Stata)

`csestudy` is available from **SSC**, and this is the recommended way to install it:

```stata
ssc install csestudy, replace
```

To update later:

```stata
ado update csestudy
```

### Bleeding-edge release from GitHub

If you want the latest release before it propagates to SSC, install directly from the GitHub Pages mirror:

```stata
net install csestudy, from("https://malcolmwardlaw.github.io/csestudy/") all replace
```

For the Python and R packages, see [python/README.md](python/README.md) and [r/README.md](r/README.md).

## Sample Data

A synthetic dataset with realistic CRSP-like properties is included for testing. It contains 300 firms over 461 S&P 500 trading days (Jan 2007 – Nov 2008), with a size-dependent abnormal return (0.15% per sd of `lag_LNMV`) injected on 2008-10-06. The panel is long enough to support GLS estimation with a 200-day pre-event window.

Load directly from GitHub:

```stata
use "https://raw.githubusercontent.com/MalcolmWardlaw/csestudy/release/examples/sample_data.dta", clear
```

Or from a local clone:

```stata
use examples/sample_data.dta, clear
```

A CSV copy (`examples/sample_data.csv`) is also included for the Python and R packages. The DGP scripts (`examples/generate_sample_data.py` and `examples/generate_sample_data.do`) are included if you want to inspect or modify the data-generating process.

## Example

```stata
* Load sample data
use "https://raw.githubusercontent.com/MalcolmWardlaw/csestudy/release/examples/sample_data.dta", clear

* Create business calendar and set panel
bcal create trading, from(date) gen(trading_date) center(20081006) replace
tsset permno trading_date

* OLS with time-series corrected errors
csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)

* GLS with 100 principal components (Cholesky, default)
csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100)

* GLS with Woodbury identity (faster, slightly less precise)
csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100) woodbury
```

## Citation

If you use this software, please cite:

```bibtex
@article{cohn2026past,
  title={Past is Prologue: Inference from the Cross Section of Returns Around an Event},
  author={Cohn, Jonathan B. and Johnson, Travis L. and Liu, Zack and Wardlaw, Malcolm I.},
  journal={Journal of Financial Economics},
  volume={180},
  pages={104278},
  year={2026},
  doi={10.1016/j.jfineco.2026.104278}
}
```

---

*Keywords: cross-sectional event study, abnormal returns, event study standard errors, time-series inference, GLS, FGLS, PCA covariance, cross-correlation in stock returns, empirical p-value, Stata, Python, R, finance, econometrics, asset pricing.*
