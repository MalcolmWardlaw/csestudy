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

{marker examples}{...}
{title:Examples}

{phang}{cmd:. bcal create trading, from(date) gen(trading_date) center( 19911121) replace}{p_end}
{phang}{cmd:. tsset permno trading_date}{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)}{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100)}{p_end}

{pstd}GLS with Woodbury identity (faster, slightly less precise):{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100) woodbury}{p_end}

{pstd}Multi-day event window using cumulative returns:{p_end}

{phang}{cmd:. gen ret5 = ret + f1.ret + f2.ret + f3.ret + f4.ret}{p_end}

{phang}{cmd:. csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-204) lastpreeventdate(-5)}{p_end}

