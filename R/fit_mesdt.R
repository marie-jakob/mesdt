#' Fit a (mixed-effects) signal detection theory (SDT) model
#'
#' @param discriminability `formula` specifying fixed and random effects on
#'  discriminability with the common syntax for mixed-effects models from
#'  `lme4` (see Details).
#' @param bias `formula` specifying fixed and random effects on response bias
#'  with the common syntax for mixed-effects models from lme4 (see Details).
#' @param dv `character` string specifying the name of the (binary) dependent
#'  variable.
#' @param trial_type `character` specifying name of the variable coding
#'  whether a trial in the given df is a signal or a noise trial.
#' @param data `data.frame` containing the variables used in the formulas for
#'  discriminability and response bias.
#' @param distribution `character` specifying the parametric distribution
#'  of signal and noise evidence ("gaussian", which is the default, "logistic",
#'  or "gumbel-min").
#' @param correlate_sdt_params `boolean` indicating whether correlations between
#'  SDT parameters should be modeled (see Details).
#' @param control list containing optional control arguments that are included
#' in the `glmmTMB()` or `glmer()` call (see Details).
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
#' @details
#' + Formulas for sensitivity and response bias can be two-sided, specifying the
#'  SDT parameter on the left-hand-side (e.g., `discriminability ~ 1 + (1 | id)`)
#'  or one-sided (e.g., `~ 1 + (1 | id)`). In the latter case, the corresponding
#'  SDT parameter is extracted from the name of the function argument.
#'
#'  + Formulas can also specify single-level models without random effects. In
#'  such a case, the `glm()` function from `stats` is used to estimate the model.
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
#'  + TODO: correlations between random effects
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
#'    bias ~ committee * emp_gender + (committee | id),
#'    data = debi3_sub,
#'    trial_type = "status",
#'    dv = "assessment"
#' )
#' summary(mod_mixed)
#'
#' # Fixed-effects model (not recommended for this type of nested data structure!):
#' mod_fixed_only <- fit_mesdt(
#'   discriminability ~ committee * emp_gender,
#'    bias ~ committee * emp_gender,
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
#' @export
fit_mesdt <- function(discriminability,
                      bias,
                      dv,
                      trial_type = "trial_type",
                      data,
                      distribution = c("gaussian", "logistic", "gumbel-min", "gumbel-max"),
                      correlate_sdt_params = T,
                      # tests = "Wald",
                      control = NULL) {
  #### Check input
  if (typeof(discriminability) != "language") stop("'discriminability' must be a formula'.")
  if (typeof(bias) != "language") stop("'bias' must be a 'formula'.")

  forms <- standardize_fit_formulas(discriminability, bias)
  discriminability <- forms[[1]]; bias <- forms[[2]]

  if (typeof(dv) != "character") stop("'dv' must be of type 'character'.")
  if (is.null(data[[dv]])) stop(paste("Given dependent variable", dv, "not in data."))
  if (length(unique(data[[dv]])) != 2) stop("dv must be a binary variable.")
  if (! all(sort(unique(data[[dv]])) == c(0, 1))) {
    if (inherits(data[[dv]], "factor") & length(unique(data[[dv]]) == 2)) {
      data[["dv_num"]] <- as.numeric(data[[dv]]) -1
      dv <- "dv_num"
    } else {
      stop("dv must be coded as 0 ('noise' response) and 1 ('signal' response)")
    }
  }
  trial_type_var <- trial_type

  if (typeof(trial_type_var) != "character") stop("'trial_type_var' must be of type 'character'.")
  if (is.null(data[[trial_type_var]])) stop(paste("Given trialtype variable", trial_type_var, "not in data."))
  # TODO: if you have a predictor that only affects sensitivity (such as strength in the context of
  # memory, this won't work) -> maybe allow a ternary variable then (maybe with a warning)
  if (all(sort(unique(data[[trial_type_var]])) != c(-1, 1))) {
    stop("'trial_type' must be a numeric binary variable coding signal trials with 1 and noise trials with -1.")
  }
  if (class(data[[trial_type_var]]) != "numeric")
    stop("'trial_type' must be of type numeric.")

  if (typeof(correlate_sdt_params) != "logical") stop("'correlate_sdt_params' must be of type 'logical'.")

  #distribution <- standardize_dist_input(distribution)
  #if (is.null(distribution)) stop("Distribution must be gaussian, logistic, gumbel-min, or gumbel-max.")
  distribution <- match.arg(distribution)

  #### Prep & fit model
  mm_all <- construct_modelmatrices(discriminability, bias, data, trial_type_var, distribution)
  m_frames <- mm_all[["m_frames"]]
  mm <- mm_all[["mm"]]

  if (distribution == "gumbel-min") {
    glmer_formula <- construct_glmer_formula(discriminability, bias, dv = "dv_rev", mm = mm,
                                             correlate_sdt_params = correlate_sdt_params)
  } else {
    glmer_formula <- construct_glmer_formula(discriminability, bias, dv = dv, mm = mm,
                                             correlate_sdt_params = correlate_sdt_params)
  }

  # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
  # (y, ID)

  fit_obj <- fit_glmm(glmer_formula, data, mm, distribution, dv, control)
  # TODO: random effects post-processing -> for the summary method

  # Check backend stuff
  if (! is.null(summary(fit_obj)$objClass[1])) {
    print("lme4 was used to fit the model.")
    backend <- "lme4"
  } else if (any(inherits(fit_obj, "glm"))) {
    print("glm() was used to fit the model.")
    backend <- "glm"
  } else {
    print("glmmTMB was used to fit the model.")
    backend <- "glmmTMB"
  }


  obj <- new_mesdt_fit(list(
    "fit_obj" = fit_obj,
    "user_input" = list(
      "discriminability" = discriminability,
      "bias" = bias,
      "dv" = dv,
      "distribution" = distribution,
      "backend" = backend,
      "trial_type_var" = trial_type_var,
      "correlate_sdt_params" = correlate_sdt_params
    ),
    "internal" = list(
      "mm" = mm,
      "m_frames" = m_frames,
      "glmer_formula" = glmer_formula,
      "data" = data,
      "backend" = backend
    )
    ))
  # Give a warning if mean sensitivity is < 0
  check_sensitivity(obj)

  return(obj)
}
