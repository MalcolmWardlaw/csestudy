# CSESTUDY: Efficient Inference for Cross-Sectional Event Studies
This is the public repository for the Stata command **csestudy** as described in Cohn, Johnson, Liu, and Wardlaw (2023) "Past is Prologue: Inference from the Cross Section of Returns Around an Event".

https://ssrn.com/abstract=4296657

The Stata program is still in beta, but it will correctly estimate the models described in the paper under reasonably general conditions and should provide a reasonable amount of error handling for the user.

Feedback is both welcome and encouraged, so please feel free to open an issue if something appears to fail or work incorrectly.

## Basic Description

As described in Cohn, Johnson, Liu, and Wardlaw (2023), testing the cross-sectional valuation effects of a specific event for firms with different characteristics is somewhat complicated. Standard event study methodologies usually fail to account for the strong cross-correlation structure in stock returns across a host of characteristics, and the standard approach of clustering the standard errors by industry is completely unable to account for this problem.

The paper proposes two additional approaches which leverage the time-series of past returns to account for the implied correlation structure:

The first is a time-series adjusted portfolio approach to inference about standard errors in which the coefficients are compared against a pre-event window of daily returns and adjusted rejection criteria are computed in the form of a parameterized z-score and a p-value estimated from the empirical distribution (the preferred metric in this approach.)

The second is a GLS based approach in which the covariance matrix Ω is also estimated from the pre-event window using a principal components approach to reduce the dimensionality. This second method adds substantial efficiency, and is the preferred approach.

## Syntax and Usage

The data must first be properly `tsset` by id and time. Further, for the default options to work, the time id must be specified as a _sequential_ integer in which non-data days like holidays and weekends are ommitted, i.e. if Friday is 10 and there are never observations on Saturday or Sunday then the following Monday is 11. The simplest way to do this is to call `bcal create` on the panel before executing the command. This method is strongly preferred as it allows the user to specify dates in a number of different ways, and the user can conveniently center the event date at t=0. See the stata help for more detail. 

The syntax is given as follows:

```stata
csestudy depvar indepvars [if], EVENTdate
  [
  EVENTENDdate(string)   ///
  NPREeventdays(integer) ///
  STARTpreeventdate(string) ENDpreeventdate(string) ///
  noBALance gls npc(integer 100) PRESAMPLEmarker(name) ///
  newvar(name) PRECALCulated
  ] 
```
The only _required_ option is:
- `EVentdate` This is the date of the event.

If the event date is specified but no other options are specified, then the pre-event window is assumed to be 200 periods long, to begin 201 periods before the event, and end 1 period before the event.

- If `NPREeventdays(n)` is specified, it overrides the length of the pre-event window default to n
- If `ENDpreeventdate()` is specified, it specifies the last date of the pre-event period, with n total pre-event days
- If `STARTpreeventdate()` is specified, it overrides the the default number of event days in favor of the first date of the pre-event period. This option cannot be specified with NPREeventdays.

- If `EVENTENDdate()` is specified, then the event is assumed to be a multi-day event and a cumulative rolling sum of `depvar` is calculated as the new dependent variable. Note that this feature is new and may still be somewhat buggy, so use with caution.
- `newvar(name)` will store the value of the cumulative rolling sum of `depvar` as a new variable specified by `name`.
- If `PRECALCulated` is specified, the program will assume that the exsting value of `depvar` is already correctly pre-calculated according to the length of the event window and will proceed as if each window is reported on the first day.
 
Additional Options:
- `npc(integer)`  Number of Principal Components. Defaults to 100
- `gls` Estimate the results using the GLS.
- `PRESAMPLEmarker(name)`  Create a variable `name` which marks the valid pre-event period observations which were used in the estimation.
- `noBALance` omits the initial balancing routine, which may slightly speed up the estimation if the user knows the sample is correctly balanced.
 

### Data Input
Data from both the event window and the pre-event window should be loaded into Stata when performing the estimation. Note that any conditional statement given by `[if]` applies *only* to the event date observations and not to any other observations.


### Balancing the pre-period data
The GLS estimation requires a strongly balanced panel in the pre-period in order to work, so any ids which do not have a full set of available returns in the pre-period will be dropped. This is done for the user by keeping only the ids which have the maximum number of observations in the pre-period. This is usually not a major issue in daily stock market data, but if your sample is massively cut down by this operation, you may have an unusual set of pre-period observations. The user should check that the data is at least *mostly* balanced before proceeding.


### Event Date Input
Note that the command will accept dates either as integer values or an Stata function which can be evaluated upon execution. Trading dates are assumed to be contiguous, but when using dates created by the `bcal` option in Stata, the command can evaluate a bcal specified date such as `eventdate(bofd("mycal",mdy(9,19,2011)))`

## Example
...
