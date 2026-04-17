/*
    generate_sample_data.do
    -----------------------
    Generates a synthetic daily stock-return panel for testing csestudy.

    DGP:
        r_{it} = alpha_i + beta_i' F_t + epsilon_{it} + abnormal_{it}

        F_t ~ N(mu_F, Sigma_F)          3 factors (market, size, value)
        beta_i ~ N(beta_bar, diag(s^2)) firm-specific loadings
        epsilon_{it} ~ N(0, sigma_i^2)  idiosyncratic, sigma_i ~ LogNormal

    Trading dates from market.stbcal (S&P 500 trading calendar).
    Event day 0 = 2008-10-06. Pre-event: 440 trading days. Post-event: 20.

    Injected event: size-dependent abnormal return on 2008-10-06.
    AR_i = 0.15% * z(lag_LNMV_i). Large firms get positive AR, small negative.

    Output:  examples/sample_data.dta
             examples/sample_data.csv

    Note: Stata's RNG differs from numpy's, so this produces a statistically
    equivalent but not numerically identical dataset to the Python version.
    The pre-generated files shipped with the repo come from the Python script.

    Usage (Stata):
        do examples/generate_sample_data.do

    Then test with:
        use examples/sample_data.dta, clear
        bcal create trading, from(date) gen(trading_date) center(20081006) replace
        tsset permno trading_date
        csestudy ret lag_LNMV if abs(prc) > 5, ///
            eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)
*/

clear all
set seed 20260417

* ── Parameters ───────────────────────────────────────────────────────────
local n_firms     = 300
local ar          = 0.0015     // abnormal return scale (per 1 sd of lag_LNMV)
local event_date  = mdy(10, 6, 2008)
local n_penny     = 15         // ~5% of firms get penny-stock prices

* ── Step 1: Build trading calendar from stbcal holidays ──────────────────
* S&P 500 trading dates, 2006-01-03 to 2009-12-31
* (Reproduces the market.stbcal omit list)

quietly {
    * Enumerate all calendar days in range
    local start_date = mdy(1, 3, 2006)
    local end_date   = mdy(12, 31, 2009)
    local total_cal  = `end_date' - `start_date' + 1

    set obs `total_cal'
    gen date = `start_date' + _n - 1
    format date %td

    * Drop weekends
    drop if dow(date) == 0 | dow(date) == 6

    * Drop market holidays
    local holidays ///
        16jan2006 20feb2006 14apr2006 29may2006 04jul2006 ///
        04sep2006 23nov2006 25dec2006  ///
        01jan2007 02jan2007 15jan2007 19feb2007 06apr2007 ///
        28may2007 04jul2007 03sep2007 22nov2007 25dec2007 ///
        01jan2008 21jan2008 18feb2008 21mar2008 26may2008 ///
        04jul2008 01sep2008 27nov2008 25dec2008 ///
        01jan2009 19jan2009 16feb2009 10apr2009 25may2009 ///
        03jul2009 07sep2009 26nov2009 25dec2009


    foreach h of local holidays {
        drop if date == td(`h')
    }

    sort date
    gen int day_index = _n

    * Find event date index and set window
    summ day_index if date == `event_date', meanonly
    local event_idx = r(mean)
    local start_idx = `event_idx' - 440
    local end_idx   = `event_idx' + 20

    keep if inrange(day_index, `start_idx', `end_idx')

    * Re-index from 1
    drop day_index
    gen int day_index = _n
    local n_days = _N

    * Verify event day
    local event_day_index = 441  // 440 pre-event + 1 (1-indexed)
    assert date == `event_date' if day_index == `event_day_index'

    tempfile calendar
    save `calendar'

}

* ── Step 2: Generate firm characteristics ────────────────────────────────

quietly {
    clear
    set obs `n_firms'
    gen int permno = 10000 + _n

    * Lagged log market value: N(8.5, 1.5) in log-millions
    * Drawn first so factor loadings can depend on size
    gen double lag_LNMV = round(rnormal(8.5, 1.5), 0.0001)

    * Factor loadings — correlated with size (realistic: small firms have
    * higher market beta, load more on size/value factors). This correlation
    * is what makes GLS outperform OLS: common factor shocks create spurious
    * day-to-day variation in the OLS coefficient on lag_LNMV that GLS removes.
    summ lag_LNMV, meanonly
    gen double lnmv_c = lag_LNMV - r(mean)
    gen double beta_market = 1.4 - 0.05 * lag_LNMV + rnormal(0, 0.20)
    gen double beta_size   = 0.6 - 0.04 * lnmv_c   + rnormal(0, 0.20)
    gen double beta_value  = 0.3 - 0.02 * lnmv_c   + rnormal(0, 0.20)
    drop lnmv_c

    * Idiosyncratic volatility (log-normal, median ~1.6%)
    gen double sigma_idio = exp(rnormal(ln(0.016), 0.35))

    * Daily alpha (small, centered near zero)
    gen double alpha = rnormal(0.0001, 0.0003)

    * Initial price: log-normal centered ~$40
    * First n_penny firms get penny-stock prices ($1.50–$4.50)
    gen double prc_initial = exp(rnormal(ln(40), 0.7))
    replace prc_initial = 1.5 + runiform() * 3 if _n <= `n_penny'

    * SIC codes: draw from non-financial industries
    local sic_list 201 283 357 366 367 382 384 481 489 737 738 739 ///
                   131 211 262 281 308 331 341 355 371 421 451 531 ///
                   541 581 701 781 871 874
    local n_sic : word count `sic_list'
    gen int sic3 = .
    forvalues i = 1/`n_firms' {
        local draw = ceil(runiform() * `n_sic')
        local code : word `draw' of `sic_list'
        replace sic3 = `code' in `i'
    }

    * Standardized lag_LNMV for size-dependent abnormal return
    summ lag_LNMV, meanonly
    local lnmv_mean = r(mean)
    summ lag_LNMV
    local lnmv_sd = r(sd)
    gen double lnmv_z = (lag_LNMV - `lnmv_mean') / `lnmv_sd'

    * Treatment indicator: above-median size
    summ lag_LNMV, detail
    gen byte treated = (lag_LNMV >= r(p50))

    tempfile firms
    save `firms'
}


* ── Step 3: Generate factor returns (time-series) ────────────────────────

quietly {
    use `calendar', clear

    * Draw 3 independent standard normals
    gen double z1 = rnormal()
    gen double z2 = rnormal()
    gen double z3 = rnormal()

    * Cholesky of correlation matrix:
    *   Corr = [[1, .2, .15], [.2, 1, .3], [.15, .3, 1]]
    *   L = [[1, 0, 0],
    *        [0.2, 0.9798, 0],
    *        [0.15, 0.2762, 0.9494]]
    local L22 = sqrt(1 - 0.2^2)
    local L32 = (0.3 - 0.2*0.15) / `L22'
    local L33 = sqrt(1 - 0.15^2 - `L32'^2)

    * Correlated standard normals → factor returns
    gen double f_market = 0.0003 + 0.011 * z1
    gen double f_size   = 0.00012 + 0.005 * (0.2 * z1 + `L22' * z2)
    gen double f_value  = 0.00012 + 0.005 * (0.15 * z1 + `L32' * z2 + `L33' * z3)

    * S&P proxy = market factor
    gen double sprtrn = round(f_market, 0.000001)

    drop z1 z2 z3

    tempfile factors
    save `factors'
}


* ── Step 4: Cross firms × dates ──────────────────────────────────────────

quietly {
    use `firms', clear
    cross using `factors'
}


* ── Step 5: Generate returns ─────────────────────────────────────────────

quietly {
    * Systematic component
    gen double systematic = beta_market * f_market ///
                          + beta_size   * f_size   ///
                          + beta_value  * f_value

    * Idiosyncratic
    gen double epsilon = rnormal() * sigma_idio

    * Total return
    gen double ret = alpha + systematic + epsilon

    * Inject size-dependent abnormal return on event day (all firms)
    * AR_i = ar * z(lag_LNMV_i) — large firms get positive, small negative
    replace ret = ret + `ar' * lnmv_z if day_index == 441

    * Round to 6 decimal places
    replace ret = round(ret, 0.000001)
}


* ── Step 6: Simulate price paths ────────────────────────────────────────

quietly {
    sort permno day_index
    by permno: gen double prc = prc_initial * exp(sum(ln(1 + ret)))
    replace prc = round(prc, 0.01)
}


* ── Step 7: Clean up and label ───────────────────────────────────────────

keep permno date ret prc sprtrn lag_LNMV sic3 treated
order permno date ret prc sprtrn lag_LNMV sic3 treated
sort permno date

label variable permno   "Firm identifier"
label variable date     "Trading date"
label variable ret      "Daily return (decimal)"
label variable prc      "Closing price (synthetic)"
label variable sprtrn   "Market index return (decimal)"
label variable lag_LNMV "Lagged log market value"
label variable sic3     "3-digit SIC code"
label variable treated  "Event treatment indicator"

label data "Synthetic CRSP-like panel for csestudy testing (seed=20260417)"

compress


* ── Step 8: Save ─────────────────────────────────────────────────────────

save examples/sample_data.dta, replace
export delimited using examples/sample_data.csv, replace


* ── Summary ──────────────────────────────────────────────────────────────

di as text _n "── Summary ────────────────────────────────────────────"
distinct permno
distinct date
count
di as text "  Event date:  " %td `event_date'
sum ret, detail
di as text "  Penny stocks (mean |prc| < $5): "
count if abs(prc) < 5
di as text "  Injected AR: +/-0.15% * z(lag_LNMV) on event day (size-dependent)"
di as text _n "  Stata usage:"
di as text "    use examples/sample_data.dta, clear"
di as text "    bcal create trading, from(date) gen(trading_date) center(20081006) replace"
di as text "    tsset permno trading_date"
di as text `"    csestudy ret lag_LNMV if abs(prc)>5, ///"'
di as text "        eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)"
