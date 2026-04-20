{smcl}
{* *! version 1.2.2  15may2018}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "csestudy##syntax"}{...}
{viewerjumpto "Description" "csestudy##description"}{...}
{viewerjumpto "Options" "csestudy##options"}{...}
{viewerjumpto "Remarks" "csestudy##remarks"}{...}
{viewerjumpto "Examples" "csestudy##examples"}{...}
{title:Title}

{phang}
{bf:csestudy} {hline 2} Efficient Inference for Cross-Sectional Event Studies


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:csestudy}
{depvar} [{indepvars}]
[{help if:if}]
{cmd:,} {opth event:startdate(csestudy##eventdate:date)} {opth firstpre:eventdate(csestudy##firstpreeventdate:date)} {opth lastpre:eventdate(csestudy##lastpreeventdate:date)} [{help csestudy##options:options}] {p_end}


{synoptset 27 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{p2coldent:* {opth event:startdate(csestudy##eventstartdate:date)}}The start date of the event.{p_end}
{p2coldent:* {opth firstpre:eventdate(csestudy##firstpreeventdate:date)}}The first (i.e. earliest) date in the pre-event period. {p_end}
{p2coldent:* {opth lastpre:eventdate(csestudy##lastpreeventdate:date)}}The last (i.e. latest) date in the pre-event period. {p_end}
{synopt:{opt gls}}Calculate GLS estimates.{p_end}
{synopt:{opt npc(integer)}}Number of principal components. Defaults to 100{p_end}
{synopt:{opt woodbury}}Use the Woodbury matrix identity for GLS instead of Cholesky decomposition. Faster but slightly less numerically precise. Requires {opt gls}.{p_end}
{synopt:{opt coefsonly}}Calculates only the coefficients, skipping the significance tests. Programmer option only.{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opth event:startdate(csestudy##eventstartdate:date)}}, {opth firstpre:eventdate(csestudy##firstpreeventdate:date)}}, and {opth lastpre:eventdate(csestudy##lastpreeventdate:date)}} are required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csestudy} calculates robust inference for cross-sectional event studies as described in {browse "https://doi.org/10.1016/j.jfineco.2026.104278":Cohn, Johnson, Liu, and Wardlaw (2026) "Past is Prologue: Inference from the Cross Section of Returns Around an Event," {it:Journal of Financial Economics} 180, 104278}.{p_end}

{pstd}
The estimation uses a time-series adjusted portfolio approach to inference about standard errors in which the coefficients are compared against a pre-event window of daily returns and adjusted rejection criteria are computed in the form of a parameterized z-score and a p-value estimated from the empirical distribution (the preferred metric in this approach.)
{p_end}

{marker options}{...}
{title:Options}


{dlgtab:Main}

{phang}
{marker eventstartdate}{...}
{opt eventstartdate(date)} The start date of the event. This date refers to the time variable set by tsset.
{p_end}

{phang}
{marker lastpreeventdate}{...}
{opt lastpreeventdate()} The last date in the pre-event period. This must be earlier than the eventstartdate later than lastpreeventdate.
{p_end}

{phang}
{marker firstpreeventdate}{...}
{opt firstpreeventdate(date)} The first date in the pre-event period. This must be earlier than the eventstartdate and lastpreeventdate.
{p_end}

{phang}
{opt gls} Calculate GLS estimates.
{p_end}

{phang}
{opt npc()} Number of principal components. Defaults to 100.
{p_end}

{phang}
{opt woodbury} Use the Woodbury matrix identity to compute the GLS transformation instead of a
Cholesky decomposition of the full covariance matrix. This inverts only a k x k matrix
(k = npc) rather than factoring the N x N covariance matrix, yielding a ~50-66% speedup
per iteration. The tradeoff is slightly reduced numerical precision due to the wide range
of idiosyncratic variances. Requires {opt gls}.
{p_end}

{marker remarks}{...}
{title:Remarks}

{pstd}
The GLS estimation requires a strongly balanced panel in the pre-period in order to work, so any ids which do not have a full set of available returns in the pre-period will be dropped. This is done for the user, and the observations which satisfy this condition are stored in e(sample). This is usually not a major issue in daily stock market data, but if your sample is significantly cut down by this operation, you may have an unusual set of pre-period observations.
{p_end}

{pstd}
Calculating significance with the estimates also require that there is a sufficiently long window of available data prior to the firstpreeventdate. (Effectively a window equal to {it:eventstartdate} - {it:firstpreeventdate} prior to firstpreeventdate). The user should check that the data is at least {it:mostly} balanced before proceeding.
{p_end}

{dlgtab:Event Date Input}

{pstd}
The date options accept integer values or Stata expressions evaluated at runtime. For example:{p_end}

{phang}{cmd:. csestudy ret lag_LNMV, eventstartdate(100) ...}{p_end}
{phang}{cmd:. csestudy ret lag_LNMV, eventstartdate(bofd("mycal",mdy(9,19,2011))) ...}{p_end}

{pstd}
It is important that the time variable tracks {it:trading days} rather than calendar days.
The simplest way to achieve this is with {help bcal:bcal create}, which constructs a
business calendar from the dates in your data and generates a sequential trading-day
variable:{p_end}

{phang}{cmd:. bcal create trading, from(date) gen(trading_date) center(20081006) replace}{p_end}
{phang}{cmd:. tsset permno trading_date}{p_end}

{pstd}
Alternatively, many users have historically used a simple sequential integer for the time
variable:{p_end}

{phang}{cmd:. bysort permno (date): gen time = _n}{p_end}
{phang}{cmd:. tsset permno time}{p_end}

{pstd}
This works as long as the panel is strongly balanced (i.e. all observations begin at the
same date and exist on all trading days). The {cmd:bcal} approach is more robust to missing
data and non-trading days, while the sequential integer approach can be more convenient if
your data is already structured that way.{p_end}

{marker sampledata}{...}
{dlgtab:Sample Data}

{pstd}
A synthetic dataset with realistic CRSP-like properties is included for testing. It contains
300 firms over 461 S&P 500 trading days (Jan 2007 {hline 1} Nov 2008), with a size-dependent
abnormal return (0.15% per sd of {bf:lag_LNMV}) injected on 2008-10-06. The panel is long
enough to support GLS estimation with a 200-day pre-event window.{p_end}

{pstd}Load directly from GitHub:{p_end}

{phang}{cmd:. use "https://raw.githubusercontent.com/MalcolmWardlaw/csestudy/release/examples/sample_data.dta", clear}{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Load the sample data and set up the panel:{p_end}

{phang}{cmd:. use "https://raw.githubusercontent.com/MalcolmWardlaw/csestudy/release/examples/sample_data.dta", clear}{p_end}
{phang}{cmd:. bcal create trading, from(date) gen(trading_date) center(20081006) replace}{p_end}
{phang}{cmd:. tsset permno trading_date}{p_end}

{pstd}OLS with time-series corrected errors:{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)}{p_end}

{pstd}GLS with 100 principal components (Cholesky, default):{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100)}{p_end}

{pstd}GLS with Woodbury identity (faster, slightly less precise):{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100) woodbury}{p_end}

{pstd}Multi-day event window using cumulative returns:{p_end}

{phang}{cmd:. gen ret5 = ret + f1.ret + f2.ret + f3.ret + f4.ret}{p_end}

{phang}{cmd:. csestudy ret5 lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-204) lastpreeventdate(-5)}{p_end}

