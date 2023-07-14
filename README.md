# CSESTUDY: Efficient Inference for Cross-Sectional Event Studies
This is the public repository for the Stata command **csestudy** as described in Cohn, Johnson, Liu, and Wardlaw (2023) "Past is Prologue: Inference from the Cross Section of Returns Around an Event".

https://ssrn.com/abstract=4296657

The Stata program is still in beta, but it will correctly estimate the models described in the paper under reasonably general conditions and should provide a reasonable amount of error handling for the user.

Feedback is both welcome and encouraged, so please feel free to open an issue if something appears to fail or work incorrectly.

## Basic Description

As described in Cohn, Johnson, Liu, and Wardlaw (2023), testing the cross-sectional valuation effects of a specific event for firms with different characteristics is somewhat complicated. Standard event study methodologies usually fail to account for the strong cross-correlation structure in stock returns across a host of characteristics, and the standard approach of clustering the standard errors by industry is completely unable to account for this problem.

The paper proposes two related approaches which leverage the time-series of past returns to account for the implied correlation structure:

The first is a time-series adjusted portfolio approach (TSOLS) in which the coefficients are compared against a pre-event window of daily returns and adjusted rejection criteria are computed in the form of an adjusted standard error, a parameterized z-score, and a p-value estimated from the empirical distribution (the preferred metric in this approach.)

The second is a GLS based approach in which the covariance matrix Ω is estimated from the pre-event window using a principal components approach to reduce the dimensionality. This second method is substantially more efficient, and is the default method provided by this command.

## Syntax and Usage

The data must first be properly `tsset` by id and time.

The syntax is given as follows:

`csestudy depvar indepvars [if] [in], EVentdate(string) STARTpreeventdate(string) ENDpreeventdate(string) [options]`

The required arguments are:
- `EVentdate` This is the date of the event.
- `STARTpreeventdate` The first date of the pre-event period
- `ENDpreeventdate` The last date of the pre-event period

Optional arguments
- `npc(integer)`  Number of Principal Components. Defaults to 100
- `tsols` Estimate the results using the TS-OLS portfolio approach to estimating the standard errors. The default is to use the GLS approach.
- `PRESAMPLEmarker(name)`  Create a variable `name` which marks the valid pre-event period observations which were used in the estimation.
- `nobalcheck` Suppresses a check to guarantee that the pre-event period is a balanced panel. This check is somewhat computationally expensive, and users may wish to supress it if they are running multiple simulations.

### Data Input
Data from both the event window and the pre-event window should be loaded into Stata when performing the estimation. The independent variables only need to be present and non-missing for the event-date observation and are completely ignored in the pre-event window. Importantly, this means that any conditional statement given by `[if]` or `[in]` applies *only* to the event date observations and not to any other observations.


### Balancing the pre-period data
The estimation requires a strongly balanced panel in the pre-period in order to work, so any ids which do not have a full set of available returns in the pre-period will be dropped. This is done for the user by keeping only the ids which have the maximum number of observations in the pre-period. This is usually not a major issue in daily stock market data, but if your sample is massively cut down by this operation, you may have an unusual set of pre-period observations. The user should check that the data is at least *mostly* balanced before proceeding.


### Event Date Input
Note that the command will accept dates either as integer values or as date functions such as `mdy()` or `td()`. Trading dates are assumed to be contiguous, but the command is fault tolerant and will accept data which has been `tsset` to include weekends and non-trading days but will also not check if unintended gaps exist where trading days did occur. 

If the cumulative returns are calculated over a window, only the date which contains full event returns should be used, and care should be taken to make sure that the earliest event window date does not overlap with the last pre-event window date.

## Example
...
