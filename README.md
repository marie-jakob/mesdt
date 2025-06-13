
<!-- README.md is generated from README.Rmd. Please edit that file -->

Estimate and test mixed-effects signal detection theory models with
maximum likelihood estimation. This package leverages the equivalence
between a subclass of SDT and a subclass of generalized linear models
first shown by DeCarlo (1998) to estimate mixed SDT models using
software for generalized linear mixed models, that is, the `lme4` or
`glmmTMB` package. Mixed-effects SDT models can be specified on the
level of SDT parameters discriminability and response bias using the
typical R formula syntax:

`discriminability = ~ X + (X | ID)`

`bias = ~ X + (X | ID)`

`mesdt` translates the given model formulas to a GLMM, uses either
`lme4` or `glmmTMB` to estimate the model (based on the user-specified
backend), and transforms the parameters back to SDT space, such that all
post-processing (e.g., estimating marginal means) can take place on the
level of SDT parameters as well. Hypothesis tests related to fixed and
random effects can be conducted with Wald tests (only for fixed
effects), likelihood ratio tests (type II and type III), and tests based
on parametric bootstrapping (type II and type III).

## Getting Started

### Step 0: Installation & Setup

You can install the development version of `mesdt` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("marie-jakob/mesdt")
```

Since `mesdt` uses `lme4` as the default backend, you also need a
working installation of `lme4` and, if you want to use it, `glmmTMB`.

``` r
library(mesdt)
```

Upon loading the package, it tells us that it has set the backend, that
is, the package used to estimate the GLMMs to `lme4`. If `glmmTMB`
(which can be significantly faster) is installed, the backend can be
changed like so:

``` r
options("mesdt.backend" = "glmmTMB")
```

### Step 1: Specify and Fit a Mixed-Effects SDT Model

We are using the `debi3` dataset provided with this package, where we
investigated attributions to gender discrimination using a signal
detection approach. In the experiment, male and female participants (as
indicated by the variable `participant gender`) had to judge 256
fictional pay raise decisions as biased or unbiased. The cases involved
male and female employees(`emp_gender`), who were either granted or
denied a pay raise (`committee`).

``` r
debi3 <- mesdt::debi3
```

To estimate effects of these variables on participants’ discriminability
and response bias, we include them as fixed effects in our model. To
account for differences between our participants in these parameters, we
include by-participant random intercepts and slopes on discriminability
and response bias. Our studies with this paradigm consistently showed
considerable by-participant variability in mean response bias, mean
discriminability and the effect of the committee decision on
participants’ response bias, which is why we include by-participant
random intercepts on both SDT parameters and a by-participant random
slope on response bias.

`mesdt` also requires the user to specify the variable indicating the
type of trial (i.e., whether the given trial was a signal, coded as 1,
or noise trial, coded as 0, `status` in `debi3`) and the variable
indicating participants’ response (`assessment` in `debi3`).

Thus, in `mesdt`, we can specify and fit our model like this:

``` r

mod <- fit_mesdt(
  discriminability ~ committee * emp_gender * participant_gender + (1 | id),
  bias ~ committee * emp_gender * participant_gender + (committee | id),
  data = debi3,
  trial_type_var = "status",
  dv = "assessment"
)
#> [1] "lme4 was used to fit the model."

summary(mod)
#> Mixed-effects signal detection theory model with Gaussian evidence distributions fit by maximum likelihood  (Adaptive Gauss-Hermite Quadrature, nAGQ = 0)with the lme4 package. 
#>  
#> Discriminability: ~committee * emp_gender * participant_gender + (1 | id) 
#> Response Bias:      ~committee * emp_gender * participant_gender + (committee | id) 
#> 
#> Fixed effects and Wald tests for discriminability: 
#>                                             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)                                 1.740768   0.051873  33.558  < 2e-16
#> committee1                                 -0.064500   0.021468  -3.004 0.002661
#> emp_gender1                                 0.029145   0.020047   1.454 0.146003
#> participant_gender1                        -0.175371   0.051873  -3.381 0.000723
#> committee1:emp_gender1                      0.003390   0.020047   0.169 0.865721
#> committee1:participant_gender1              0.013531   0.021468   0.630 0.528510
#> emp_gender1:participant_gender1             0.035473   0.020047   1.769 0.076812
#> committee1:emp_gender1:participant_gender1 -0.009961   0.020047  -0.497 0.619261
#> 
#> Fixed effects and Wald tests for response bias: 
#>                                             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)                                 0.151158  -0.026793   5.642 1.68e-08
#> committee1                                  0.031319  -0.076754   0.408    0.683
#> emp_gender1                                -0.015726  -0.010023  -1.569    0.117
#> participant_gender1                        -0.023609  -0.026793  -0.881    0.378
#> committee1:emp_gender1                      0.016032  -0.010024   1.599    0.110
#> committee1:participant_gender1             -0.005700  -0.076754  -0.074    0.941
#> emp_gender1:participant_gender1            -0.001055  -0.010023  -0.105    0.916
#> committee1:emp_gender1:participant_gender1  0.011400  -0.010024   1.137    0.255
```

The `summary()` method prints population-level estimates for the fixed
effects we specified separate for discriminability and response bias in
a similar format as`lme4` and `glmmTMB`. The p values in the last column
are results from Wald tests based on the beta estimates and their
standard errors.

### Step 2: Compute Hypothesis Tests for Selected Parameters

Since Wald tests are often not recommended for inference in GLMM,
`mesdt` allows the user to compute different types (see the
documentation for details) of likelihood ratio tests and parametric
bootstrapping tests, that are based on comparisons of nested models.

`mesdt` allows the user to either specify the to-be-tested parameters
via formulas (e.g., `~ committee`) or to test all fixed effects
(`"all"`). For this experiment, we were interested (among other things)
if participants’ response bias and / or discriminability vary as a
function of the committee decision:

``` r
fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                 trial_type_var = "trial_type_fac")
#> [1] "lme4 was used to fit the model."

lrts_test <- compute_tests(fit, test_intercepts = T)

tests <- compute_tests(mod, 
                       tests = "lrt",
                       tests_discriminability = ~ committee,
                       tests_response_bias = ~ committee)

tests$LRTs$LRT_results
#>                  deviance_full deviance_reduced df.LRT Chisq     p.value    
#> committee_lambda 21584.14      21584.31         1      0.1608541 0.6883714  
#> committee_mu     21584.14      21593.07         1      8.924664  0.002813451
```

In this subset of our data, there is a significant effect of the
committee decision on discriminability, but not on response bias.

### Step 3: Estimate marginal means

`mesdt` additionally provides a custom `emmeans()` function, allowing
the user to estimate marginal means for response bias and
discriminability. The syntax is the typical `emmeans` syntax with an
additional argument `dpar` specifying the SDT parameter. Thus, we can
get estimated marginal means for response bias and discriminability for
“denied” and “granted” decisions like so:

``` r
library(emmeans)

emmeans(mod, ~ committee, dpar = "discriminability")
#> NOTE: Results may be misleading due to involvement in interactions
#>  committee emmean     SE  df asymp.LCL asymp.UCL
#>  true        1.68 0.0558 Inf      1.57      1.79
#>  false       1.81 0.0565 Inf      1.69      1.92
#> 
#> Results are averaged over the levels of: emp_gender, participant_gender 
#> Confidence level used: 0.95
emmeans(mod, ~ committee, dpar = "response bias")
#> NOTE: Results may be misleading due to involvement in interactions
#>  committee emmean     SE  df asymp.LCL asymp.UCL
#>  true       0.182 0.0785 Inf    0.0286     0.336
#>  false      0.120 0.0840 Inf   -0.0448     0.284
#> 
#> Results are averaged over the levels of: emp_gender, participant_gender 
#> Confidence level used: 0.95
```

The estimated marginal means show that participants’ discriminability
was higher in trials where a pay raise was denied than for cases where
one was granted. Descriptively,

## References

TODO
