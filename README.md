
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

Upon loading the package, it tells us that it has set the backend, that
is, the package used to estimate the GLMMS to `lme4`. If `glmmTMB`
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

To estimate effects of these variables on participants’ sensitivity and
response bias, we include them as fixed effects in our model. To account
for differences between our participants in these parameters, we include
by-participant random intercepts and slopes on sensitivity and response
bias. Our studies with this paradigm consistently showed considerable
by-participant variability in mean response bias, mean sensitivity and
the effect of the committee decision on participants’ response bias,
which is why we include by-participant random intercepts on both SDT
parameters and a by-participant random slope on response bias.

`mesdt` also requires the user to specify the variable indicating the
type of trial (i.e., whether the given trial was a signal, coded as 1,
or noise trial, coded as 0, `status` in `debi3`) and the variable
indicating participants’ response (`assessment` in `debi3`).

Thus, in `mesdt`, we can specify and fit our model like this:

``` r

# TODO: change input so that lhs formulas are accepted as well 
# (kind of unintuitive like this)
mod <- fit_mesdt(
  discriminability =~ committee * emp_gender * participant_gender + (1 | id),
  bias =~ committee * emp_gender * participant_gender + (committee | id),
  data = debi3,
  trial_type_var = "status",
  dv = "assessment"
)
#> [1] "glmmTMB was used to fit the model."
```

### Step 2: Compute Hypothesis Tests for Selected Parameters

### Step 3: Estimate marginal means

## References

- glmmTMB package
- lme4 package
- DeCarlo 1998
