
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
the number of interim analyses (**nia**) and potentially time point
(**tia**). the distribution (**simfun**) and paremeters in each arm
(**par01**), the estimation function (**estfun**),

Let’s assume

- a trial with a target final sample size of 200 (1:1 allocation),
- a ‘good’ binary outcome, which is expected to occur in 20% of the
  control patient,
- a relevant effect would be an absolute increase by 20%,
- a risk difference as effect measure using Wald-type CI and test.

We would like to do one interim analysis at half of patients Single
simulation under the null and alternative can the be done like that:

``` r
set.seed(12)

simH0<-ADsim(n01 = c(100, 100), nia = 1, 
    par01 = c(p0=0.2, p1=0.2), simfun = "binom",
    estfun = "rd_wald", 
    direct = "higher")

simH1<-ADsim(n01 = c(100, 100), nia = 1, 
    par01 = c(p0=0.2, p1=0.4), simfun = "binom",
    estfun = "rd_wald", 
    direct = "higher")
    
cbind(simH0,simH1)
#>                simH0        simH1
#> true_pe   0.00000000   0.20000000
#> fa_n0   100.00000000 100.00000000
#> fa_n1   100.00000000 100.00000000
#> fa_x0    17.00000000  21.00000000
#> fa_x1    23.00000000  42.00000000
#> fa_p0     0.17000000   0.21000000
#> fa_p1     0.23000000   0.42000000
#> fa_pe     0.06000000   0.21000000
#> fa_lci   -0.05056004   0.08457762
#> fa_uci    0.17056004   0.33542238
#> fa_z      1.06365592   3.28165062
#> ia1_n0   50.00000000  50.00000000
#> ia1_n1   50.00000000  50.00000000
#> ia1_x0   10.00000000   7.00000000
#> ia1_x1   15.00000000  24.00000000
#> ia1_p0    0.20000000   0.14000000
#> ia1_p1    0.30000000   0.48000000
#> ia1_pe    0.10000000   0.34000000
#> ia1_lci  -0.06860248   0.17139752
#> ia1_uci   0.26860248   0.50860248
#> ia1_z     1.16247639   3.95241972
```

`ADsim` returns a numeric vector with the true point estimate
(*true_pe*, which is a risk difference, a log risk ratio or a mean
difference for the default options), and—for each interim stage (prefix
*iax\_*) and the final stage (prefix *fa\_*)—the

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
    ADsim(n01 = c(100, 100), nia = 1, 
    par01 = c(p0=0.2, p1=0.2), simfun = "binom",
    estfun = "rd_wald", 
    direct = "higher"))
simdatH0 <- do.call(rbind, simdatH0)

simdatH1 <- lapply(1:1000, function(i) 
    ADsim(n01 = c(100, 100), nia = 1, 
    par01 = c(p0=0.2, p1=0.4), simfun = "binom",
    estfun = "rd_wald", 
    direct = "higher"))
simdatH1 <- do.call(rbind, simdatH1)
```

### User-defined estimation function

Estimation functions can be specified via **estfun**. There are five
defaults, standard Wald-type approximation for the risk difference
(*wald_rd*) and risk ratio (*wald_rr*), the score-based method for the
risk difference suggested by Newcombe (*score_rd*), the score-based
method for the risk ratio according to Koopman (*score_rr*), and
Student’s t-test for the mean difference (*md_t*).

Alternatively, any function can be defined, which takes hte output from
**simfun** and *cilevel* as arguements.

The output of the function must be a (named) vector with point estimate,
lower confidence limit, upper confidence limit and z-statistic. It can
also contain further elements (e.g. estimates for each group), which are
added to the output of `ADsim`.

If a user-specified function is given in **estfun**, option **truefun**
has also to be used. It must be a function with argument *par01* that
specifies how the true effect is calculated, e.g. as *diff(par01)* for
the risk difference.

This example would reproduce the *wald_rd* option.

``` r

efun<-function(out,cilevel) {
    xs<-unlist(lapply(out,sum))
    ns<-unlist(lapply(out,length))
    pe<-diff(xs/ns)
    se<-sqrt(sum(xs/ns*(1-xs/ns)/ns))
    lci<-pe - qnorm(1-(1-cilevel)/2)*se
    uci<-pe + qnorm(1-(1-cilevel)/2)*se
    z<-(pe)/se
    res<-c(pe,lci,uci,z)
    names(res)<-c("pe","lci","uci","z")
    return(res)
}

tfun <- function(x) diff(x)

ADsim(n01 = c(100,100), par01 = c(0.2,0.4), nia = 1, estfun = efun, truefun = tfun)
```

### Different outcome types

Continuous outcomes can be specified like that:

``` r

par01<-data.frame(matrix(c(5,1,4.8,1),2,2))
rownames(par01)<-c("mean","sd")
ADsim(n01 = c(100,100), par01 = par01, nia = 1, simfun="norm", estfun = "md_t")
```

Note that the *par01* needs rows named *mean* and *sd*.

Other outcome types can also be included using specific functions for
**simfun**, **estfun** and **truefun**.

**simfun** has to take par01 as argument, and produce a list of length
two with the outcomes for control and experimental group.

For example, a Poisson with expected counts of 4 and 3 in control and
experimental group and low counts beeing better:

``` r

sfun<-function(n01,par01) {
    out0<-rpois(n = n01[1], lambda = par01[1])
    out1<-rpois(n = n01[2], lambda = par01[2])
    return(list(out0,out1))
}

efun<-function(out,cilevel) {
    
    ns<-unlist(lapply(out,length))
    xs<-unlist(lapply(out,mean))
    
    m1<-glm(c(out[[1]],out[[2]]) ~ c(rep(0,ns[1]),rep(1,ns[2])),poisson(link = "log"))
    
    z<-summary(m1)$coef[2,"z value"]
    
    res<-c(xs,m1$coef[2],confint.default(m1, level = cilevel)[2,],z)

    names(res)<-c("m0","m1","pe","lci","uci","z")
    return(res)
}

tfun<-function(p01) diff(log(par01))

ADsim(n01 = c(100,100), nia = 1,
    par01 = c(4,3), simfun=sfun,
    estfun = efun, truefun = tfun,
    direct="lower")
```

## Evaluation

Repeated simulations can be evaluated using function `ADeval`.

Thresholds for stopping at each interim analysis have to be defined via
**zcutoff**, a matrix with one row per interim analysis and a column for
stopping for futility and efficacy. A stop for futility will be assumed
if the observed Z-statistic at interim is lower than the threshold in
column 1, a stop for efficacy if the observed Z-statistic at interim is
higher than the threshold in column 2. *NA* would mean no stopping at
the respective interim analysis.

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

head(resH0$rawdata)
#>   true_pe fa_n0 fa_n1 fa_x0 fa_x1 fa_p0 fa_p1 fa_pe      fa_lci     fa_uci       fa_z ia1_n0 ia1_n1 ia1_x0 ia1_x1 ia1_p0 ia1_p1 ia1_pe     ia1_lci    ia1_uci      ia1_z obs_n0 obs_n1 obs_x0 obs_x1
#> 1       0   100   100    17    23  0.17  0.23  0.06 -0.05056004 0.17056004  1.0636559     50     50     10     15   0.20   0.30   0.10 -0.06860248 0.26860248  1.1624764    100    100     17     23
#> 2       0   100   100    21    18  0.21  0.18 -0.03 -0.13974048 0.07974048 -0.5357997     50     50      7     10   0.14   0.20   0.06 -0.08677500 0.20677500  0.8012116    100    100     21     18
#> 3       0   100   100    20    17  0.20  0.17 -0.03 -0.13754828 0.07754828 -0.5467212     50     50     12      6   0.24   0.12  -0.12 -0.26875081 0.02875081 -1.5811388     50     50     12      6
#> 4       0   100   100    17    19  0.17  0.19  0.02 -0.08645329 0.12645329  0.3682298     50     50      6      7   0.12   0.14   0.02 -0.11177023 0.15177023  0.2974821    100    100     17     19
#> 5       0   100   100    12    24  0.12  0.24  0.12  0.01481730 0.22518270  2.2360680     50     50      9     10   0.18   0.20   0.02 -0.13372916 0.17372916  0.2549892    100    100     12     24
#> 6       0   100   100    14    16  0.14  0.16  0.02 -0.07893451 0.11893451  0.3962144     50     50      5     11   0.10   0.22   0.12 -0.02176922 0.26176922  1.6590038    100    100     14     16
#>   obs_p0 obs_p1 obs_pe     obs_lci    obs_uci      obs_z tstop tstop_fut tstop_eff anystop_fut anystop_eff              anystop
#> 1   0.17   0.23   0.06 -0.05056004 0.17056004  1.0636559    NA        NA        NA       FALSE       FALSE          Not stopped
#> 2   0.21   0.18  -0.03 -0.13974048 0.07974048 -0.5357997    NA        NA        NA       FALSE       FALSE          Not stopped
#> 3   0.24   0.12  -0.12 -0.26875081 0.02875081 -1.5811388     1         1        NA        TRUE       FALSE Stopped for futility
#> 4   0.17   0.19   0.02 -0.08645329 0.12645329  0.3682298    NA        NA        NA       FALSE       FALSE          Not stopped
#> 5   0.12   0.24   0.12  0.01481730 0.22518270  2.2360680    NA        NA        NA       FALSE       FALSE          Not stopped
#> 6   0.14   0.16   0.02 -0.07893451 0.11893451  0.3962144    NA        NA        NA       FALSE       FALSE          Not stopped

cbind(resH0$opchar,resH1$opchar)
#>                       [,1]         [,2]
#> true_pe         0.00000000   0.20000000
#> av_n0          75.45000000  68.40000000
#> av_n1          75.45000000  68.40000000
#> av_x0          15.07100000  13.65300000
#> av_x1          15.03400000  27.48900000
#> av_p0           0.20640000   0.19256000
#> av_p1           0.19273000   0.41194000
#> av_pe          -0.01367000   0.21938000
#> av_lci         -0.14553059   0.06507062
#> av_uci          0.11819059   0.37368938
#> av_z           -0.10639400   2.76331804
#> psig            0.04500000   0.89400000
#> pstop           0.49100000   0.63200000
#> pstop_fut       0.45900000   0.01300000
#> pstop_eff       0.03200000   0.61900000
#> pstop_ia1       0.49100000   0.63200000
#> pstop_fut_ia1   0.45900000   0.01300000
#> pstop_eff_ia1   0.03200000   0.61900000
#> av_n          150.90000000 136.80000000
#> min_n         100.00000000 100.00000000
#> max_n         200.00000000 200.00000000
#> bias           -0.01367000   0.01938000
#> sd              0.07469009   0.08056865
#> mse             0.00575990   0.00686040
#> coverage        0.91500000   0.92900000
```

`ADeval` returns a list with the evaluations for each repetition
(*rawdata*) and the summarized operating characteristics (*opchar*).

*rawdata* is a data frame with the elements from `ADsim` plus

- prefix *obs\_*: the observed sample size (*n0*, *n1*), effect (*pe*)
  and confidence interval (*lci*, *uci*) at study end (which could be at
  an interim or the final stage), and potentially further variables
  included in `ADsim::estfun`.
- *tstop*: the time point of the first stop for either futility or
  efficacy,  
  and of the first stop for futility and efficacy (*tstop_fut* and
  *tstop_eff*).
- *anystop*: an indicator whether there was any stop, and a stop for
  futility or efficacy (*anystop_fut* and *anystop_eff*).

*opchar* is a numeric vector with operating characteristics over all
simulated trials including the

- prefix *av\_*: average estimates in each group, effect (point
  estimate, *pe*) and 95% confidence limits (*lci*, *uci*).
- *psig*: probability for a significant trial (the power or type I
  error, depending on the assumptions).
- *pstop*: probability for any stop and a stop for futility or efficacy
  (*pstop_fut*, *pstop_eff*), overall and at each interim analysis
  (suffix *iax*).
- average, minimum and maximum total sample size (av_n, min_n, max_n).
- *bias*: average of estimate minus true value).
- *sd*: standard deviation of the point estimate.
- *mse*: mean squared error, average of the squared difference between
  estimate and true value).
- *coverage*: proportion of confidence intervals including the true
  value.

## shinyApp

A shinyApp can be launched locally in RStudio via `launch_ADSim_app()`.
