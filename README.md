
**Beta Software Disclaimer**  
\> This package is in beta. Functionality may change, and some features
are experimental.  
\> Use at your own discretion and report issues on the GitHub Issues
page.

# mesdt

<!-- README.md is generated from README.Rmd. Please edit that file -->

This package leverages the equivalence between a subclass of SDT and a
subclass of generalized linear models first shown by DeCarlo (1998) to
estimate SDT models using software for (generalized) linear mixed models
(GLMM), that is, the `lme4` or `glmmTMB` package (for mixed-effects SDT
models) and the `glm()` function (for single-level SDT models). SDT
models can be specified on the level of the SDT parameters
discriminability and response bias using the typical R formula syntax:

`discriminability = ~ X + (X | ID)`

`response_bias = ~ X + (X | ID)`

`mesdt` translates the given model formulas to a generalized linear
(mixed) model, uses either `lme4`, `glmmTMB` or `glm()` to estimate the
model (based on the user-specified backend), and transforms the
parameters back to SDT space, such that all post-processing (e.g.,
estimating marginal means) can take place on the level of SDT parameters
as well. Hypothesis tests can be conducted with Wald tests (returned be
the fitting function for each fixed effects parameter), likelihood ratio
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
working installation of `lme4` and, if you want to use it, `glmmTMB`. We
are also loading the `tibble` package for a nicer display of the data.

``` r
library(mesdt)
# install.packages("tibble")
library(tibble)
```

Upon loading the package, it tells us that it has set the backend, that
is, the package used to estimate the GLMMs, to `lme4`. If `glmmTMB`
(which can be significantly faster) is installed, the backend can be
changed like so:

``` r
set_backend("glmmTMB")
```

### Step 1: Prepare the Data

We are using the `debi3subset` dataset provided with this package[^1],
where we investigated attributions to gender discrimination using a
signal detection approach. In the experiment, participants had to judge
256 fictional pay raise decisions as biased or unbiased. The cases
involved male and female employees(`emp_gender`), who were either
granted or denied a pay raise (`committee`).

``` r
debi3subset <- as_tibble(mesdt::debi3subset)
```

`mesdt` assumes the data to be in the long format, meaning that one row
of the data represents one observation of the binary response variable.

Fitting an SDT model requires at least two variables: the binary
response variable (`assessment` in `debi3subset`) and the type of trial
(i.e., whether the given trial was a signal or a noise trial, `status`
in `debi3subset`). Both can be factors or numeric variables, where 1 (or
in case of a factor, the level coded as one in the contrast coding)
corresponds to a signal and 0 or -1 corresponds to a noise response or
trial.

``` r
head(debi3subset)
#> # A tibble: 6 × 10
#>   id    assessment status objective committee emp_gender file_name      participant_gender   age status_num
#>   <fct> <fct>      <fct>  <chr>     <fct>     <fct>      <fct>          <fct>              <dbl>      <dbl>
#> 1 1     unbiased   noise  true      granted   m          W_M_20.png     m                     35         -1
#> 2 1     unbiased   noise  true      granted   m          W_M_47.png     m                     35         -1
#> 3 1     unbiased   noise  true      granted   m          W_M_f_4.jpg    m                     35         -1
#> 4 1     biased     signal true      denied    f          W_F_o_u_10.jpg m                     35          1
#> 5 1     biased     signal true      denied    f          W_F_122.png    m                     35          1
#> 6 1     biased     signal true      denied    m          W_M_10.png     m                     35          1
```

For categorical predictor variables, we strongly recommend using sum
contrasts for a more straightforward interpretation of lower-order
effects in the presence of higher-order terms. We also strongly
recommend centering continuous variables.

``` r
contrasts(debi3subset$committee) <- contr.sum(2)
contrasts(debi3subset$emp_gender) <- contr.sum(2)
```

### Step 2: Specify and Fit a Mixed-Effects SDT Model

To estimate effects of the employee’s gender and the decision of the
committee on participants’ discriminability and response bias, we
include them as fixed effects in our model. To account for differences
between our participants in these parameters, we include by-participant
random intercepts and slopes on discriminability and response bias. For
the present paradigm, the maximal random-effects structure justified by
the design includes both by-participant random intercepts and random
slopes for the effects of the committee decision, employee gender and
their interaction on both discriminability and response bias. To
considerable reduce the complexity of our model and speed up estimation,
we do not estimate correlations between the random effects here, by
dividing the random effects and the grouping factor with two vertical
bars (`||`). Note that unlike lme4 and glmmTMB, mesdt also supports this
notation for factors.

Thus, in `mesdt`, we can specify and fit our model like this:

``` r

mod <- fit_mesdt(
  discriminability =~ committee * emp_gender + (committee * emp_gender || id),
  response_bias =~ committee * emp_gender + (committee * emp_gender || id),
  data = debi3subset,
  trial_type = "status",
  dv = "assessment"
)
#> Given random-effects structure contains uncorrelated terms. Modeling all random effects parametes as uncorrelated since a mix of correlated and uncorrelated terms is not supported at the moment.
#> Correlating SDT Parameters is not possible in the presence of uncorrelated terms.
#> boundary (singular) fit: see help('isSingular')

summary(mod)
#> Mixed-effects signal detection theory model with Gaussian evidence distributions fit by maximum likelihood (Laplace Approximation) with the lme4 package. 
#>  
#> Discriminability:  ~committee * emp_gender + (committee * emp_gender || id) 
#> Response Bias:     ~committee * emp_gender + (committee * emp_gender || id) 
#> Data:  debi3subset 
#> 
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>    4753.8    4858.5   -2360.9    4721.8      5104 
#> 
#> Random effects:
#>  Groups Name                                     Std.Dev. 
#>  id     (Intercept)(Response Bias)               3.876e-01
#>  id.1   committee1(Response Bias)                5.718e-01
#>  id.2   emp_gender1(Response Bias)               0.000e+00
#>  id.3   committee1:emp_gender1(Response Bias)    1.076e-05
#>  id.4   (Intercept)(Discriminability)            3.853e-01
#>  id.5   committee1(Discriminability)             0.000e+00
#>  id.6   emp_gender1(Discriminability)            8.369e-02
#>  id.7   committee1:emp_gender1(Discriminability) 0.000e+00
#> Number of obs: 5120, groups:  id, 20
#> 
#> Fixed effects and Wald tests for discriminability: 
#>                        Estimate Std. Error z value Pr(>|z|)
#> (Intercept)             1.77862    0.09986  17.811   <2e-16
#> committee1              0.04763    0.04706   1.012    0.311
#> emp_gender1            -0.01445    0.04729  -0.306    0.760
#> committee1:emp_gender1 -0.03823    0.04335  -0.882    0.378
#> 
#> Fixed effects and Wald tests for response bias: 
#>                         Estimate Std. Error z value Pr(>|z|)
#> (Intercept)             0.141268   0.089565   1.577    0.115
#> committee1              0.015486   0.130010   0.119    0.905
#> emp_gender1             0.005583   0.021799   0.256    0.798
#> committee1:emp_gender1 -0.035406   0.021855  -1.620    0.105
#> 
#> optimizer (Nelder_Mead) convergence code: 0 (OK)
#> boundary (singular) fit: see help('isSingular')
```

In practice, we recommend to iteratively reduce the random-effects
structure such that only random-effects terms that are supported by the
data are included in the model (see Bates et al., 2015 and Jakob et al.,
2025 for explanations).

The `summary()` method prints population-level estimates for the fixed
effects and variances and covariances of the random effects separately
for discriminability and response bias in a similar format as `lme4` and
`glmmTMB`. The p values in the last column are results from Wald tests
based on the beta estimates and their standard errors.

### Step 3: Compute Hypothesis Tests for Selected Parameters

Since Wald tests are often not recommended for inference in GLMM,
`mesdt` allows the user to compute different types (see the
documentation for details) of likelihood ratio tests and parametric
bootstrapping tests, that are based on comparisons of nested models via
the `compute_tests()` function.

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
#> Correlating SDT Parameters is not possible in the presence of uncorrelated terms.
#> Correlating SDT Parameters is not possible in the presence of uncorrelated terms.
#> boundary (singular) fit: see help('isSingular')
#> boundary (singular) fit: see help('isSingular')

tests
#> Type III likelihood ratio tests 
#> 
#> Discriminability: 
#>           deviance_full deviance_reduced df.LRT Chisq   p.value
#> committee       4721.82          4722.83      1  1.01 0.3152676
#> 
#> Response Bias: 
#>           deviance_full deviance_reduced df.LRT Chisq   p.value
#> committee       4721.82          4721.83      1  0.01 0.9045006
```

In this subset of our data, there is neither a significant effect of the
committee decision on discriminability, nor on response bias.

### Step 4: Estimate marginal means

`mesdt` additionally provides a custom `emmeans()` function, allowing
the user to estimate marginal means for response bias and
discriminability. The syntax is the typical
[`emmeans`](https://cran.r-project.org/web/packages/emmeans/index.html)
syntax with an additional argument `dpar` specifying the SDT parameter.
Thus, we can get estimated marginal means for response bias and
discriminability for “denied” and “granted” decisions like so:

``` r
library(emmeans)

emmeans(mod, ~ committee, dpar = "discriminability")
#> NOTE: Results may be misleading due to involvement in interactions
#>  committee emmean    SE  df asymp.LCL asymp.UCL
#>  denied      1.83 0.110 Inf      1.61      2.04
#>  granted     1.73 0.111 Inf      1.51      1.95
#> 
#> Results are averaged over the levels of: emp_gender 
#> Confidence level used: 0.95
emmeans(mod, ~ committee, dpar = "response bias")
#> NOTE: Results may be misleading due to involvement in interactions
#>  committee emmean    SE  df asymp.LCL asymp.UCL
#>  denied     0.157 0.158 Inf    -0.153     0.466
#>  granted    0.126 0.158 Inf    -0.184     0.435
#> 
#> Results are averaged over the levels of: emp_gender 
#> Confidence level used: 0.95
```

The estimated marginal means show that participants’ discriminability
was descriptively lower in trials where a pay raise was denied
($d' = 1.73$) than for cases where one was granted ($d' = 1.83$). We can
see a different pattern for response bias, with a smaller, that is, more
liberal response criterion in “granted” cases ($c = 0.126$) than in
“denied” cases ($c = 0.157$). However, both effects were not
significantly different from zero in this subset of the data.

[^1]: `debi3subset` contains a subset of 20 participants. The complete
    dataset is provided as `debi3`.
