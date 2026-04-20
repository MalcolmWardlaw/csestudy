/*
    tests/test_help_examples.do
    ---------------------------
    Executes the Examples section from csestudy.sthlp to verify the
    command runs end-to-end against the shipped sample data.

    Run from the repo root:
        do tests/test_help_examples.do

    `set varabbrev off` is set per the SSC archive submission guidelines
    (http://repec.org/bocode/s/sscsubmit.html): programs must work
    regardless of the user's varabbrev preference.
*/

clear all
set more off
set varabbrev off

* Use the repo's ado files rather than any SSC-installed copy
adopath ++ "`c(pwd)'"

* ── Load sample data ────────────────────────────────────────────────────
* (The help file loads from a GitHub URL; the local copy is used here so
*  the test is self-contained.)
use examples/sample_data.dta, clear

* Create business calendar and set panel
bcal create trading, from(date) gen(trading_date) center(20081006) replace
tsset permno trading_date

* ── Example 1: OLS with time-series corrected errors ────────────────────
csestudy ret lag_LNMV if abs(prc)>5, ///
    eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)

* ── Example 2: GLS with 100 principal components (Cholesky, default) ────
csestudy ret lag_LNMV if abs(prc)>5, ///
    eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) ///
    gls npc(100)

* ── Example 3: GLS with Woodbury identity ───────────────────────────────
csestudy ret lag_LNMV if abs(prc)>5, ///
    eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) ///
    gls npc(100) woodbury

* ── Example 4: Multi-day event window using cumulative returns ──────────
gen ret5 = ret + f1.ret + f2.ret + f3.ret + f4.ret
csestudy ret lag_LNMV if abs(prc)>5, ///
    eventstartdate(0) firstpreeventdate(-204) lastpreeventdate(-5)

di as result _n "All help-file examples executed successfully."
