# Change-point detection and homogenisation of daily climatic time series.

Implementation of the Standard Normal Homogeneity Test (SNHT, Alexandersson 1986) to detect and correct change points in daily time series of environmental variables.

## Citation

If you use this or part of this code in your work, please cite it as:

Beguería, S. 2023. Change-point detection and homogeneisation of daily time series of climatic data, version 1.0.0. https://github.com/sbegueria/homogen [date of last access].

## Details

There are three main functions available (file 'functions_snht.R'):

* `snht()` detects a single change point in a test series and tests for its statistical significance.
* `snht.Q()` can be used to build a test series suited to apply the SNHT from one candidate time series and several auxiliary ones.
* `snht.daily()` uses the two functions above to apply the SNHT to a daily time series. The function applies the procedure recursively and corrects the found inhomogeneities, until no more significant change points are found.

The first two functions leverage on code by 

There are two auxiliary functions:

* `snht.plot()` produces a diagnostic plot with the output of a call to `snht()`.
* `snht.daily.plot()` produces diagnostic plots with the output of a call to `snht.daily()`.

In the notebook we apply the procedure to time series of minimum temperature.


## Version history

### Version 1.0.0, April 2023.

First version on GitHub.
