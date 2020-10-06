ad hoc methods for detecting population trends
================
David Nguyen
2020-10-06

Trend detection: sequential vs fixed sample approaches
======================================================

Wildlife managers and conservationists need to detect population trends in a timely and reliable manner. Timeliness is reflected in the number of samples (e.g., years) needed to detect a trend. Reliability is captured by the type I error and power of a test. Simple linear regression using sampling time as a covariate is commonly used to detect population trends. Researchers have used this approach to identify the minimum years of monitoring needed to detect population trends.

For example, [White (2019)](https://academic.oup.com/bioscience/article/69/1/40/5195956#129750432) repeatedly fit linear models to sub-samples of a time series to identify the minimum sample size needed to achieve a specified type I error and power. That is, for a time series with *T* time points, he fit a linear regression model to each contiguous sub-samples of size 2, 3, …, *T* − 1. Using the arbitrary (but conventional) benchmark of 0.8 power at the 0.05 significance level, White used the regression outputs to find the minimum sub-sample length such that 80 % of the samples had a significant regression slope coefficient. This is essentially a sample size calculation for a one-sided hypothesis test with *H*<sub>0</sub> : *r* = 0 and *H*<sub>1</sub> : *r* &gt; 0 or *H*<sub>1</sub> : *r* &lt; 0.

White used data simulated from the following population model for examples of how to calculate fixed-sample sizes needed to detect population trends,
*N*<sub>*t* + 1</sub> ∼ *N*(*N*<sub>*t*</sub> + *r*, *σ*),

Where *N*<sub>*t*</sub> is the population size at time *t*, *r* is the population trend, and *σ* is the population variability.

We can also test this one-sided hypothesis using a sequential test. We can define our sequential test statistic at time *t* (*S*<sub>*t*</sub>) as:

$\\hat{r}\_t$ is the MLE of r computed from data points *X*<sub>1</sub>, …, *X*<sub>*t*</sub>. The MLE is constrained to (0, inf) or ( − inf, 0) depending on whether detection of a linear increase or decrease is desired. (In my computer implementation, I just use a finite interval for the constrained MLE).

The decision rule is:

The thresholds are set following Wald's method such that $A \\sim \\frac{1-\\beta}{\\alpha}$ and $B \\sim \\frac{\\beta}{1-\\alpha}$ where *α* and *β* is the probabilities of type I and II error, respectively (I'm not sure if these approximations work for this inference problem, I still need to check if these boundaries give the correct error probabilities).

For the situation where *σ* is unknown, the sequential test statistic is instead calculated using the current ML estimate of of the standard deviation ($\\hat{\\sigma\_t}$). $\\hat{\\sigma\_t}$ is computed as a [running variance](https://www.johndcook.com/blog/standard_deviation/) from the first differences of the observations $N\_t - N\_{t-1} \\overset{iid}{\\sim} Norm(r,\\sigma)$ for *t* = 2, 3, 4, …. If we know that there is a lower bound on the variability of the population when it is stable we can use *σ*<sub>*t*</sub> = *m**a**x*(*σ*<sub>*m**l**e*</sub>,*σ*<sub>*m**i**n*</sub>).

Comparison of detection times
=============================

Tests:

1.  fixed-sample size linear regression test
2.  sequential test (*σ* known)
3.  sequential test (*σ* unknown)
    -   the running MLE estimate of *σ* is used to calculate the test statistic
    -   note that the MLE estimate of *σ* is uncorrected (biased low)

I compare the performance of the tests according to the time required to reject the null hypothesis that there is no population trend (*α* = 0.05, *β* = 0.2). The population data were simulated from an additive model *N*<sub>*t* + 1</sub> = *N*(*N*<sub>*t*</sub> + *r*, *σ*) with parameters *r*= 1.5 and *σ*= 5. The population was initialized as 1000 individuals and was run for 100 time steps. 100 simulations were used.

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)

The fixed-sample size needs 20 samples to achieve 0.8 power at the 0.05 significance level. In comparison, the sequential tests require 12 (*σ* known) and 11 (*σ* unknown and *σ*<sub>*m**i**n*</sub>= 3). The reason that the sequential test where *σ* is unknown detects the change faster than when *σ* is known (for the parameters used here) is because the MLE estimate of *σ* is uncorrected which biases it low. This makes the likelihood of the null smaller than it should be. The table belows shows the number of instances either sequential method required as many or more samples than the fixed-sample size test.

| method       | delay &gt;= samplesize |    n|
|:-------------|:-----------------------|----:|
| delay\_est   | FALSE                  |   96|
| delay\_est   | TRUE                   |    4|
| delay\_known | FALSE                  |   92|
| delay\_known | TRUE                   |    8|

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Note, the probability models used for the fixed-sample and sequential tests are not the same. The sequential test uses the true data-generating model (variability is from process noise) whereas the fixed-sample test fits a simple linear regression using time as a covariate (observation noise). A more comparable test would be to use a sequential t-test on the first-differenced values of the population sizes.
