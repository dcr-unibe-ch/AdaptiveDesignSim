
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AdaptiveDesignSim

<!-- badges: start -->

<!-- [![R-CMD-check](https://github.com/dcr-unibe-ch/AdaptiveDesignSim/workflows/R-CMD-check/badge.svg)](https://github.com/dcr-unibe-ch/AdaptiveDesignSim/actions) -->

<!-- badges: end -->

An R package for simulations of clinical trials with adaptive designs.

## Installation

The package can be installed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dcr-unibe-ch/AdaptiveDesignSim")
```

## Simulation

A single trial with one or more interim analyses can be simulated using
`ADsim`. The required options are the sample size in each arm (**n01**),
the proportion of the outcome in each arm (**p01**), and the number of
interim analyses (**nia**).

Further optional arguments for `ADsim` are

- the time point of the interim analysis (**tia**, as fraction of the
  total, equally spaced by default),
- the estimation function (**estfun**, either a user-specified function
  or a string with one of the four defaults (wald_rd, wald_rr, score_rd,
  score_rr),
- the two-sided confidence level for the confidence interval if one of
  the defaults was used for *estfun*,
- a function to calculate the true point estimate if a user-specified
  *estfun* was used (**truefun**),  
- the direction of the effect (**direct** with *lower* is better or
  *higher* is better).

Let’s assume

- a trial with a target final sample size of 200 (1:1 allocation),
- a ‘bad’ binary outcome, which is expected to occur in 40% of the
  control patient,
- a relevant effect would be a decrease by 20%,
- a risk difference as effect measure using Wald-type CI and test with a
  one-sided alpha of 0.025.

We would like to do one interim analysis at half of patients Single
simulation under the null and alternative can the be done like that:

``` r
set.seed(12)

simH0<-ADsim(n01 = c(100, 100), p01 = c(0.4, p1=0.4), nia = 1, direct = "lower")

simH1<-ADsim(n01 = c(100, 100), p01 = c(0.4, p1=0.2), nia = 1, direct = "lower")

cbind(simH0,simH1)
#>                simH0         simH1
#> true_p0   0.40000000   0.400000000
#> true_p1   0.40000000   0.200000000
#> true_pe   0.00000000  -0.200000000
#> fa_x0    36.00000000  41.000000000
#> fa_x1    43.00000000  18.000000000
#> fa_n0   100.00000000 100.000000000
#> fa_n1   100.00000000 100.000000000
#> fa_p0     0.36000000   0.410000000
#> fa_p1     0.43000000   0.180000000
#> fa_pe     0.07000000  -0.230000000
#> fa_lci   -0.06515227  -0.352321225
#> fa_uci    0.20515227  -0.107678775
#> fa_z     -1.01513261   3.685310674
#> ia1_x0   18.00000000  19.000000000
#> ia1_x1   27.00000000  10.000000000
#> ia1_n0   50.00000000  50.000000000
#> ia1_n1   50.00000000  50.000000000
#> ia1_p0    0.36000000   0.380000000
#> ia1_p1    0.54000000   0.200000000
#> ia1_pe    0.18000000  -0.180000000
#> ia1_lci  -0.01179627  -0.354337667
#> ia1_uci   0.37179627  -0.005662333
#> ia1_z    -1.83941802   2.023621877
```

`ADsim` returns a numeric vector with

- the true proportions in both arms (*true_p0*, *true_p1*).
- the true point estimate (*true_pe*), which is either a risk difference
  or a log risk ratio (depending on option **estfun**).

And for each interim stage (prefix *iax\_*) and the final stage (prefix
*fa\_*) the

- simulated number of successes for both arms (suffix *x0*, *x1*),
- sample size for both arms (suffix *n0*, *n1*),
- simulated probabilities for both arms (suffix *p0*, *p1*),
- simulated point estimate (suffix *pe*, either the risk difference or
  the log risk ratio),
- two-sided *cilevel*% confidence intervals (suffix *lci* and *uci*),
  and
- z-statistic (suffix *z*).

Simulations can be repeated and stacked to generate a data frame that
can then be used with `ADeval`.

``` r
set.seed(12)

simdatH0 <- lapply(1:1000, function(i) 
    ADsim(n01 = c(100,100), p01 = c(0.2, p1=0.2), nia = 1, 
          direct="lower"))
simdatH0 <- do.call(rbind, simdatH0)

simdatH1 <- lapply(1:1000, function(i) 
    ADsim(n01 = c(100,100), p01 = c(0.4, p1=0.2), nia = 1,
          direct="lower"))
simdatH1 <- do.call(rbind, simdatH1)
```

Optionally a different effect estimate or estimation function can be
specified via **estfun**. There are four defaults, standard Wald-type
approximation for the risk difference (*wald_rd*) and risk ratio
(*wald_rr*), the score-based method for the risk difference suggested by
Newcombe (*score_rd*) and the score-based method for the risk ratio
according to Koopman (*score_rr*).

Alternatively, any function can be defined with arguments

- *xs*: the number of successes in control and experimental group)
- *ns*: the number of patients in control and experimental group
- *cilevel*: the confidence level

and a named vector as output with point estimate, lower confidence
limit, upper confidence limit and z-statistic.

In that case, function **truefun** with argument *p01* has to be used to
specify how the true effect should be calculated, e.g. as *diff(p01)*
for the risk difference.

This example would reproduce the **wald_rd** option:

``` r

efun<-function(xs,ns,cilevel) {
    pe<-diff(xs/ns)
    se<-sqrt(sum(xs/ns*(1-xs/ns)/ns))
    lci<-pe - qnorm(1-(1-cilevel)/2)*se
    uci<-pe + qnorm(1-(1-cilevel)/2)*se
    z<-(-pe)/se
    res<-c(pe,lci,uci,z)
    names(res)<-c("pe","lci","uci","z")
    return(res)
}

ADsim(n01 = c(100,100), p01 = c(0.4,0.2), nia = 1, estfun = efun, truefun = function(x) diff(x))
```

## Evaluation

For the `ADeval` thresholds for stopping at each interim analysis have
to be defined via **zcutoff**, a matrix with one row per interim
analysis and a column for stopping for futility and efficacy. A stop for
futility will be assumed if the observed Z-statistic at interim is lower
than the threshold in column 1, a stop for efficacy if the observed
Z-statistic at interim is higher than the threshold in column 2. *NA*
would mean no stopping at the respective interim analysis.

Let’s assume we have one interim analysis (as specified in the `ADsim`
above) and would like to stop for futility if Z\<0 and for efficacy if
Z\>2:

``` r
zcutoff<-matrix(c(0,2),1,2)
zcutoff
#>      [,1] [,2]
#> [1,]    0    2
```

``` r
resH0<-ADeval(simdat = simdatH0, zcutoff = zcutoff)

resH1<-ADeval(simdat = simdatH1, zcutoff = zcutoff)

cbind(resH0$opchar,resH1$opchar)
#>                       [,1]         [,2]
#> true_p0         0.20000000   0.40000000
#> true_p1         0.20000000   0.20000000
#> true_pe         0.00000000  -0.20000000
#> av_p0           0.19317000   0.40656000
#> av_p1           0.20625000   0.19501000
#> av_pe           0.01308000  -0.21155000
#> av_lci         -0.11811740  -0.36451660
#> av_uci          0.14427740  -0.05858340
#> psig            0.04300000   0.88800000
#> pstop           0.47500000   0.60300000
#> pstop_fut       0.44800000   0.01400000
#> pstop_eff       0.02700000   0.58900000
#> pstop_ia1       0.47500000   0.60300000
#> pstop_fut_ia1   0.44800000   0.01400000
#> pstop_eff_ia1   0.02700000   0.58900000
#> avn           152.50000000 139.70000000
#> minn          100.00000000 100.00000000
#> maxn          200.00000000 200.00000000
#> bias            0.01308000  -0.01155000
#> sd              0.07287952   0.07665359
#> mse             0.00547720   0.00600330
#> coverage        0.91900000   0.93900000
```

`ADeval` returns a list with the evaluations for each repetition
(*rawdata*) and the summarized operating characteristics (*opchar*).

*rawdata* is a data frame with the elements from `ADsim` plus

- the observed proportions, effect and confidence interval at study end,
  which could be at an interim or the final stage (*obs_p0*, *obs_p1*,
  *obs_pe*, *obs_lci*, *obs_uci*),
- the time point of the first stop for either futility or efficacy
  (*tstop*),  
  and of the first stop for futility and efficacy (*tstop_fut* and
  *tstop_eff*),
- the final sample size (*neff*), and
- an indicator whether there was any stop (*anystop*) and a stop for
  futility or efficacy (*anystop_fut* and *anystop_eff*).

*opchar* is a numeric vector with operating characteristics over all
simulated trials including the

- average estimates in each group, effect (point estimate, *pe*) and 95%
  confidence limits (*lci*, *uci*),
- probability for a significant trial (*psig*, the power or type I
  error, depending on the assumptions),
- probability for any stop (*pstop*) and a stop for futility or efficacy
  (*pstop_fut*, *pstop_eff*), overall and at each interim analysis
  (suffix *iax*),
- average, minimum and maximum sample size (avn, minn, maxn),
- bias (average of estimate minus true value),
- standard deviation (*sd*) of the point estimate,
- mean squared error (*mse*, average of the squared difference between
  estimate and true value), and
- coverage (proportion of confidence intervals including the true
  value).

## shinyApp

A shinyApp can be launched locally in RStudio via `launch_ADSim_app()`.
