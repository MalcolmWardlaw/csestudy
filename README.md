# CSESTUDY: Efficient Inference for Cross-Sectional Event Studies
This is the public repository for the Stata command **csestudy** as described in Cohn, Johnson, Liu, and Wardlaw (2025) "Past is Prologue: Inference from the Cross Section of Returns Around an Event".

https://ssrn.com/abstract=4296657

The Stata program is still in very early beta, but it will correctly estimate the models described in the paper under reasonably general conditions and should provide a reasonable amount of error handling for the user.

Feedback is both welcome and encouraged, so please feel free to open an issue if something appears to fail or work incorrectly.

## Basic Description

As described in Cohn, Johnson, Liu, and Wardlaw (2025), testing the cross-sectional valuation effects of a specific event for firms with different characteristics is somewhat complicated. Standard event study methodologies usually fail to account for the strong cross-correlation structure in stock returns across a host of characteristics, and the standard approach of clustering the standard errors by industry is completely unable to account for this problem.

The paper proposes an approach which leverage the time-series of past returns to account for the implied correlation structure:

The estimation uses a time-series adjusted portfolio approach to inference about standard errors in which the coefficients are compared against a pre-event window of daily returns and adjusted rejection criteria are computed in the form of a parameterized z-score and a p-value estimated from the empirical distribution (the preferred metric in this approach.)


## Syntax and Usage

The data must first be properly **tsset** by id and time. Further, for the default options to work, the time id must be specified as a _sequential_ integer in which non-data days like holidays and weekends are ommitted, i.e. if Friday is 10 and there are never observations on Saturday or Sunday then the following Monday is 11. The simplest way to do this is to call **bcal create** on the panel before executing the command. This method is strongly preferred as it allows the user to specify dates in a number of different ways, and the user can conveniently center the event date at t=0. See the stata help for more detail. 

The syntax is given as follows:

```stata
csestudy depvar [indepvars] [if] , eventstartdate() firstpreeventdate() lastpreeventdate() 
```
*Additional Options*

```stata
gls
npc(integer)
coefsonly
```

### Data Input
Data from both the event window and the pre-event window should be loaded into Stata when performing the estimation. Note that the conditional statement given by **[if]** applies to the event date and pre-event-date observations, but not to the y variables in used for calculating the PCA matrix if the gls option is specified.


### GLS and balancing the pre-period data
The GLS estimation requires a strongly balanced panel of nonmissing values for the independent variable in each pre-period. It also requires that there is a sufficiently long window of available data prior to the firstpreeventdate. (Effectively a window equal to _eventstartdate_ - _firstpreeventdate_ prior to _firstpreeventdate_.)  The user should check that the data is at least *mostly* balanced across this window before proceeding.


### Event Date Input
Note that the command will accept dates either as integer values or an Stata function which can be evaluated upon execution. Trading dates are assumed to be contiguous, but when using dates created by the **bcal** option in Stata, the command can evaluate a bcal specified date such as `eventdate(bofd("mycal",mdy(9,19,2011)))`

## Installation
```stata
net install csestudy, from("https://malcolmwardlaw.github.io/csestudy/") all replace
```

## Example

```stata
bcal create trading, from(date) gen(trading_date) center(20081006) replace
tsset permno trading_date

csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1)
csestudy ret lag_LNMV if abs(prc)>5, eventstartdate(0) firstpreeventdate(-200) lastpreeventdate(-1) gls npc(100)

```
