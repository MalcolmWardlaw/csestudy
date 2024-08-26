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
{cmd:,} {opth event:date(csestudy##eventdate:eventdate)} [{help csestudy##options:options}] {p_end}


{synoptset 27 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{p2coldent:* {opth event:date(csestudy##eventdate:eventdate)}}The start date of the event{p_end}
{synopt:{opt eventend:date(string)}}Last date of the event. Defaults to the same day as eventdate.{p_end}
{synopt:{opt npre:eventdays(integer)}}Number of pre-event days in window. Defaults to 200.{p_end}
{synopt:{opt end:preeventdate(string)}}Last pre-event date. Defaults to one period prior to the event date.{p_end}
{synopt:{opt start:preeventdate(string)}}Start of pre-event date. If un-specified, defaults to n event days prior to the preeventdate.{p_end}
{synopt:{opt gls}}Calculate GLS estimates.{p_end}
{synopt:{opt npc(integer)}}Number of principle components. Defaults to 100{p_end}
{synopt:{opt presample:marker(newvar)}}Create variable newvar which marks the pre-event sample.{p_end}
{synopt:{opt newvar(varname)}}Preserves the multi-day calculated returns in a new variable.{p_end}
{synopt:{opt precalc:ulated}}If the multi-day returns in the depvar have been pre-caulcuated for each cell, this option will take them as given rather than calculating them on the fly.{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opth event:date(csestudy##eventdate:eventdate)}} is required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csestudy} calculates robust inference for cross-sectional event studies as described in {browse "https://ssrn.com/abstract=4296657":Cohn, Johnson, Liu, and Wardlaw (2024) "Past is Prologue: Inference from the Cross Section of Returns Around an Event"}.{p_end}

{pstd}
The estimation uses a time-series adjusted portfolio approach to inference about standard errors in which the coefficients are compared against a pre-event window of daily returns and adjusted rejection criteria are computed in the form of a parameterized z-score and a p-value estimated from the empirical distribution (the preferred metric in this approach.){p_end}

{marker options}{...}
{title:Options}


{dlgtab:Main}

{phang}
{marker eventdate}{...}
{opt eventdate(eventdate)} The start date of the event is a required option. If the event is a single period long, this is the only necessary option. The program will use the existing defaults to calculate the significance statistics using the default 200 period pre-event window. If the event is longer than a day, then {opt eventenddate()} must be specified. If it is, the program will calculate a rolling n period window of returns, making sure to skip n period before the end of the pre-event window.

{phang}
{opt eventenddate()} Last date of the event. Defaults to the same day as eventdate.

{phang}
{opt npreeventdays()} Number of pre-event days in window. Defaults to 200.

{phang}
{opt endpreeventdate()} Last pre-event date. Defaults to one period prior to the event date.

{phang}
{opt startpreeventdate()} Start of pre-event date. If un-specified, defaults to n event days prior to the preeventdate.

{phang}
{opt gls} Calculate GLS estimates.

{phang}
{opt npc()} Number of principle components. Defaults to 100

{phang}
{opt presample:marker(newvar)} Create variable newvar which marks the pre-event sample.

{phang}
{opt newvar()} Preserves the multi-day calculated returns in a new variable.

{phang}
{opt precalculated} If the multi-day returns in the depvar have been pre-caulcuated for each cell, this option will take them as given rather than calculating them on the fly.

{marker remarks}{...}
{title:Remarks}

{pstd}
The GLS estimation requires a strongly balanced panel in the pre-period in order to work, so any ids which do not have a full set of available returns in the pre-period will be dropped. This is done for the user by keeping only the ids which have the maximum number of observations in the pre-period. This is usually not a major issue in daily stock market data, but if your sample is massively cut down by this operation, you may have an unusual set of pre-period observations. The user should check that the data is at least {it:mostly} balanced before proceeding.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. bcal create trading, from(date) gen(trading_date) center(20081006) replace}{p_end}
{phang}{cmd:. tsset permno trading_date}{p_end}

{phang}{cmd:. csestudy ret btm me, eventdate(0)}{p_end}

{phang}{cmd:. csestudy ret btm me, start(-210) end(-11) eventdate(0)}{p_end}

