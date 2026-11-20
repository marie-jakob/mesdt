#' Fit a (mixed-effects) signal detection theory (SDT) model
#'
#' @param discriminability `formula` specifying fixed and random effects on
#'  discriminability with the common syntax for mixed-effects models from
#'  `lme4` (see Details).
#' @param response_bias `formula` specifying fixed and random effects on response bias
#'  with the common syntax for mixed-effects models from lme4 (see Details).
#' @param dv `character` string specifying the name of the (binary) dependent
#'  variable. Can be numeric (1 = signal response, 0 or -1 = noise response) or
#'  a factor. For factors, the level coded as 1 in the contrasts
#'  (`contrasts(data[[dv]])`) is treated as signal response, unless
#'  `dv_signal_level` is given.
#' @param trial_type `character` specifying name of the variable coding
#'  whether a trial in the given df is a signal or a noise trial. Can be numeric
#'  (1 = signal trial, 0 or -1 = noise trial) or a factor. For factors, the level
#'  coded as 1 in the contrasts (`contrasts(data[[trial_type]])`) is treated as
#'  signal trial, unless `trial_type_signal_level` is given.
#' @param data `data.frame` containing the variables used in the formulas for
#'  discriminability and response bias.
#' @param distribution `character` specifying the parametric distribution
#'  of signal and noise evidence ("gaussian", which is the default, "logistic",
#'  or "gumbel-min").
#' @param correlate_sdt_params `boolean` indicating whether correlations between
#'  SDT parameters should be modeled (default is TRUE, see Details).
#' @param aggregate (experimental) `boolean` indicating whether the given long
#'  data frame should be aggregated. Can significantly speed up estimation,
#'  especially for lme4, but has the disadvantage that additional methods based
#'  on the fitted model (e.g., bootstrapping) might not work as expected.
#' @param control list containing optional control arguments that are included
#' in the `glmmTMB()` or `glmer()` call (see Details).
#' @param na.action `character` a function that indicates what should happen
#' when the data contain NAs. The default action (na.omit, inherited from the
#' ‘factory fresh’ value of getOption("na.action")) strips any observations
#' with any missing values in any variables.
#' @param dv_signal_level optional `character` string giving the level of `dv`
#'  that codes signal responses (only if `dv` is a factor). If `NULL` (default),
#'  the level coded as 1 in the contrasts of `dv` is used.
#' @param trial_type_signal_level optional `character` string giving the level
#'  of `trial_type` that codes signal trials (only if `trial_type` is a factor).
#'  If `NULL` (default), the level coded as 1 in the contrasts of `trial_type`
#'  is used.
#'
#' @return An object of class `mesdt_fit`, containing the fitted model
#' (`$fit_obj`), information about the specified strucure of the model
#' (`$user_input`) and internal information used for post-processing.
#' (`$internal`).
#'
#' @description
#' Estimates (mixed-effects) signal detection theory (SDT) models with maximum
#' likelihood estimation by leveraging the equivalence between certain SDT and
#' certain generalized linear models (GLMM; De Carlo, 1998). The GLMMs are
#' estimated using either the `lme4` (the default) or the `glmmTMB` package
#' (which can be significantly faster) as a backend, which can be changed with
#' \link{set_backend}.
#'
#' The resulting parameter estimates can be interpreted on the latent evidence
#' scale that is assumed in SDT: For instance, when using sum contrasts, the fixed (or population-level)
#' intercepts for discriminability (d') and response bias (c) indicate mean overall
#' levels of discriminability and response bias. The higher d', the larger the
#' difference between the means of the signal and noise distributions and the
#' better signal and noise can be discriminated. A response bias of c < 0 (c > 0)
#' indicates a liberal (conservative) response tendency, that is, a tendency
#' to give a "signal" ("noise") response. The fixed effects of predictors
#' (again, assuming sum contrasts for all predictors) describe the population
#' effect of an predictor on discriminability or response bias. (i.e., for a
#' categorical variable with two levels, the difference between those levels).
#' The random effects indicate the variability between the levels of the
#' random-effects grouping factor (e.g., between participants).
#' For a more detailed explanation and tutorial see Jakob et al. (2025).
#'
#' The default SDT parametrization is the common equal-variance Gaussian model,
#' but `fit_mesdt()` and this package support other distributions (i.e., the
#' logistic, gumbel-min and gumbel-max distribution) as well, which can be
#' more appropriate in certain contexts.
#'
#' `summary()` and `print()` methods are provided, showing parameter estimates,
#'  standard errors and results of Wald tests (as returned by `lme4` and
#'  `glmmTMB`); however, it is usually recommended to rely on likelihood ratio
#' tests or parametric bootstrapping tests (in case of few levels of the
#' random-effects grouping factor) for statistical inference with GLMM, which
#' are provided in the \link{compute_tests} function.
#'
#' In addition, the returned `mesdt_fit` object contains the fitted GLMM,
#' (`$fit_obj`, a `glmerfit`), allowing users to do additional processing that is not (yet)
#' implemented within this package, such as post-processing of random effect
#' conditional modes.
#'
#' @details
#' + Formulas for discriminability and response bias should be one-sided (e.g.,
#' `discriminability =~ 1 + (1 | id)`).
#'
#'  + Suppressing correlations between random effects is possible by dividing
#'  the random-effects formula and the grouping factor with two, instead of one
#'  vertical bars (`||`). Note that unlike lme4 and glmmTMB, mesdt also supports
#'  this notation for factors.
#'
#'  + The current implementation supports completely correlated and
#'  completely uncorrelated random-effects structures, as well as correlations
#'  only within random effects for the same SDT parameter (i.e., random effects
#'  for discriminability are correlated with each other but not with random
#'  effects for response bias and vice versa). Formulas containing both
#'  correlated and uncorrelated terms are not yet supported and will be estimated
#'  as completely uncorrelated models.
#'
#'  + Formulas can also specify single-level models without random effects. In
#'  such a case, the `glm()` function from `stats` is used to estimate the model.
#'
#'
#'  + Per default, correlations between the random effects for discriminability
#'  and response bias are modeled with a joint covariance matrix (i.e., random
#'  intercepts and slopes for discriminability and response bias are correlated
#'  with each other). Setting `correlate_sdt_params = F` splits this up into
#'  two separate covariance matrices for discriminability and response bias
#'  (which can help to achieve convergence).
#'
#'  + The `control` argument allows to pass all control arguments taken by
#'  `lme4` or `glmmTMB` (depending on the backend used; see the relevant
#'  documentation for details).
#'
#'
#' @examples
#' \dontrun{
#' # Mixed-effects SDT model
#' # by-participant random intercepts for msensitivity and response bias
#' # by-participant random slope for the effect of the committee decision on
#' # response bias
#' # correlations between all random effects
#' mod_mixed <- fit_mesdt(
#'   discriminability ~ committee * emp_gender + (1 | id),
#'    response_bias ~ committee * emp_gender + (committee | id),
#'    data = debi3_sub,
#'    trial_type = "status",
#'    dv = "assessment"
#' )
#' summary(mod_mixed)
#'
#' # Fixed-effects model (not recommended for this type of multi-level data structure!):
#' mod_fixed_only <- fit_mesdt(
#'   discriminability ~ committee * emp_gender,
#'   response_bias ~ committee * emp_gender,
#'    data = debi3_sub,
#'    trial_type = "status",
#'    dv = "assessment"
#' )
#' summary(mod_fixed_only)
#' }
#' @importFrom lme4 fixef
#' @importFrom lme4 ranef
#' @importFrom lme4 VarCorr
#' @importFrom stats vcov
#' @importFrom stats contrasts
#' @export
fit_mesdt <- function(discriminability,
                      response_bias,
                      dv,
                      trial_type,
                      data,
                      distribution = c("gaussian", "logistic", "gumbel-min", "gumbel-max"),
                      correlate_sdt_params = TRUE,
                      aggregate = FALSE,
                      control = NULL,
                      na.action = c("omit", "fail"),
                      dv_signal_level = NULL,
                      trial_type_signal_level = NULL) {
  bias <- response_bias
  #### Check input
  if (typeof(discriminability) != "language") stop("'discriminability' must be a formula'.")
  if (typeof(bias) != "language") stop("'response_bias' must be a 'formula'.")

  check_fit_formulas(discriminability, bias)

  check_random_intercepts(discriminability, "discriminability")
  check_random_intercepts(bias, "response bias")


  data_name <- deparse(substitute(data))

  dv_input <- dv
  trial_type_input <- trial_type
  # Check and code dv as 0/1 (column dv_num) and trial type as -1/1 (column trial_type_num);
  # missing values stay NA and are handled below
  dv_coded <- check_dv(data, dv, dv_signal_level)
  data <- dv_coded$data
  dv <- dv_coded$var

  trial_type_coded <- check_trial_type_var(data, trial_type, trial_type_signal_level)
  data <- trial_type_coded$data
  trial_type_var <- trial_type_coded$var

  # Remove trials with missing values in any variable used by the model, so that all models
  # (including submodels in compute_tests()) are fitted to the same trials
  data <- remove_incomplete_trials(data, c(dv_input, trial_type_input,
                                           all.vars(discriminability), all.vars(bias)))

  data_input <- data


  if (typeof(correlate_sdt_params) != "logical") stop("'correlate_sdt_params' must be of type 'logical'.")

  distribution <- match.arg(distribution)

  # Check datatypes of predictors
  check_predictors(data, discriminability, response_bias)

  # reverse-code dv for gumbel-min distribution
  if (distribution == "gumbel-min") {
    data[["dv_rev"]] <- ifelse(data[[dv]] == 0, 1, 0)
  }

  if (aggregate) {
    data_orig <- data
    if (distribution != "gumbel-min") {
      data <- aggregate_data(data, discriminability, response_bias, "dv_num", trial_type_var)
    } else {
      data <- aggregate_data(data, discriminability, response_bias, "dv_rev", trial_type_var)
    }
  }

  dv <- ifelse(aggregate, "dv_agg",
               ifelse(distribution == "gumbel-min", "dv_rev", dv))


  #### Prep & fit model
  mm_all <- construct_modelmatrices(discriminability, bias, data, trial_type_var, distribution)
  m_frames <- mm_all[["m_frames"]]
  mm <- mm_all[["mm"]]

  mm <- mm_all[["mm"]]
  check_identifiable(mm)

  check_random_structure(discriminability, bias, data,
                         correlate_requested = !missing(correlate_sdt_params) && correlate_sdt_params)

  glmer_formula <- construct_glmer_formula(discriminability, bias, dv = dv, mm = mm,
                                           correlate_sdt_params = correlate_sdt_params)

  # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
  # (y, ID)

  fit_obj <- fit_glmm(glmer_formula, data, mm, distribution, dv, control)
  # TODO: random effects post-processing -> for the summary method

  # Check backend stuff
  if (! is.null(summary(fit_obj)$objClass[1])) {
    #cat("Model was estimated with lme4.")
    backend <- "lme4"
  } else if (any(inherits(fit_obj, "glm"))) {
    #cat("Model was estimated with glm().")
    backend <- "glm"
  } else {
    #cat("Model was estimated with glmmTMB.")
    backend <- "glmmTMB"
  }


  obj <- new_mesdt_fit(list(
    "fit_obj" = fit_obj,
    "user_input" = list(
      "discriminability" = discriminability,
      "bias" = bias,
      "dv_input" = dv_input,
      "data_name" = data_name,
      "data_input" = data_input,
      "distribution" = distribution,
      "backend" = backend,
      "trial_type_input" = trial_type_input,
      "dv_signal_level" = dv_coded$signal_level,
      "trial_type_signal_level" = trial_type_coded$signal_level,
      "correlate_sdt_params" = correlate_sdt_params,
      "aggregate" = aggregate
    ),
    "internal" = list(
      "mm" = mm,
      "m_frames" = m_frames,
      "glmer_formula" = glmer_formula,
      "data_mod" = data,
      "backend" = backend,
      "dv" = dv,
      "trial_type" = trial_type
    )
    ))
  # Give a warning if mean sensitivity is < 0
  check_sensitivity(obj)

  return(obj)
}
