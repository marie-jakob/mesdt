# mesdt
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

`mesdt` translates the given model formulas to a generalized linear
mixed model (GLMM), uses either `lme4` or `glmmTMB` to estimate the
model (based on the user-specified backend), and transforms the
parameters back to the SDT logic, such that all post-processing (e.g.,
estimating marginal means) can take place on the level of SDT parameters
as well. 

Hypothesis tests can be conducted with Wald tests (returned be the fitting
function for each fixed effects parameter), likelihood ratio
tests (type II and type III), and tests based on parametric
bootstrapping (type II and type III).

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
is, the package used to estimate the GLMMs, to `lme4`. If `glmmTMB`
(which can be significantly faster) is installed, the backend can be
changed like so:

``` r
set_backend("glmmTMB")
```

### Step 1: Specify and Fit a Mixed-Effects SDT Model

We are using the `debi3_sub` dataset provided with this package[^1],
where we investigated attributions to gender discrimination using a
signal detection approach. In the experiment, participants had to judge
256 fictional pay raise decisions as biased or unbiased. The cases
involved male and female employees(`emp_gender`), who were either
granted or denied a pay raise (`committee`).

``` r
debi3_sub <- mesdt::debi3_sub
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
or noise trial, coded as 0, `status` in `debi3_sub`) and the variable
indicating participants’ response (`assessment` in `debi3_sub`).

Thus, in `mesdt`, we can specify and fit our model like this:

``` r

mod <- fit_mesdt(
  discriminability ~ committee * emp_gender + (1 | id),
  bias ~ committee * emp_gender + (committee | id),
  data = debi3_sub,
  trial_type = "status",
  dv = "assessment"
)
#> [1] "lme4 was used to fit the model."

summary(mod)
#> Mixed-effects signal detection theory model with Gaussian evidence distributions fit by maximum likelihood (Adaptive Gauss-Hermite Quadrature, nAGQ = 0) with the lme4 package. 
#>  
#> Discriminability:  ~committee * emp_gender + (1 | id) 
#> Response Bias:     ~committee * emp_gender + (committee | id) 
#> 
#> Random effects:
#>  Groups Name                          Std.Dev. Corr       
#>  id     (Intercept)(Response Bias)    0.3848              
#>         committee1(Response Bias)     0.5801    0.17      
#>         (Intercept)(Discriminability) 0.3902   -0.39  0.18
#> 
#> Fixed effects and Wald tests for discriminability: 
#>                        Estimate Std. Error z value Pr(>|z|)
#> (Intercept)             1.75932    0.09974  17.638   <2e-16
#> committee1             -0.04001    0.04671  -0.857    0.392
#> emp_gender1            -0.01419    0.04308  -0.329    0.742
#> committee1:emp_gender1  0.03843    0.04308   0.892    0.372
#> 
#> Fixed effects and Wald tests for response bias: 
#>                         Estimate Std. Error z value Pr(>|z|)
#> (Intercept)             0.145662  -0.088899   1.639    0.101
#> committee1             -0.018567  -0.131821  -0.141    0.888
#> emp_gender1             0.006281  -0.021540   0.292    0.771
#> committee1:emp_gender1  0.034437  -0.021540   1.599    0.110
```

The `summary()` method prints population-level estimates for the fixed
effects and variances and covariances of the random effects separately
for discriminability and response bias in a similar format as `lme4` and
`glmmTMB`. The p values in the last column are results from Wald tests
based on the beta estimates and their standard errors.

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
tests <- compute_tests(mod, 
                       tests = "lrt",
                       tests_discriminability = ~ committee,
                       tests_response_bias = ~ committee)
#> [1] "committee_lambda" "committee_mu"

tests
#> Type III likelihood ratio tests 
#> 
#> Discriminability: 
#>           deviance_full deviance_reduced df.LRT Chisq   p.value
#> committee       4718.78          4719.49      1  0.71 0.3999445
#> 
#> Response Bias: 
#>           deviance_full deviance_reduced df.LRT Chisq   p.value
#> committee       4718.78           4718.8      1  0.02 0.8803103
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
#>  committee emmean   SE  df asymp.LCL asymp.UCL
#>  true        1.72 0.11 Inf      1.50      1.94
#>  false       1.80 0.11 Inf      1.58      2.01
#> 
#> Results are averaged over the levels of: emp_gender 
#> Confidence level used: 0.95
emmeans(mod, ~ committee, dpar = "response bias")
#> NOTE: Results may be misleading due to involvement in interactions
#>  committee emmean    SE  df asymp.LCL asymp.UCL
#>  true       0.127 0.170 Inf    -0.207     0.461
#>  false      0.164 0.147 Inf    -0.123     0.452
#> 
#> Results are averaged over the levels of: emp_gender 
#> Confidence level used: 0.95
```

The estimated marginal means show that participants’ discriminability
was higher in trials where a pay raise was denied ($d' = 1.72$ than for
cases where one was granted $d' = 1.83$). Descriptively, we can see a
similar pattern for response bias, with a smaller, that is, more liberal
response criterion in “denied” cases ($\lambda = 0.12$) than in
“granted” cases ($\lambda = 0.183$), but this difference was not
significantly different from zero in this subset of the data.

## References

TODO

[^1]: `debi3_sub` contains a subset of 20 participants. The complete
    dataset is provided as `debi3`.
