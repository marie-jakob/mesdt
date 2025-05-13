#' Fit a multilevel signal detection theory model
#'
#' @param discriminability Formula specifying fixed and random effects on sensitivity
#' @param bias Formula specifying fixed and random effects on response bias
#' @param dv name of the (binary) dependent variable
#' @param trial_type_var name of the variable coding signal vs. noise trials
#' @param data dataset
#' @param correlate_sdt_params boolean indicating whether correlations between
#'  SDT parameters should be modeled
#' @param tests type of statistical tests that should be used ("Wald": Wald tests,
#'  "LRT": likelihood ratio tests, "boot": parametric bootstrapping)
#'
#' @return TODO
#' @importFrom lme4 fixef
#' @importFrom lme4 ranef
#' @importFrom lme4 VarCorr
#' @importFrom stats vcov
#' @export
#'
#' @examples
fit_mesdt <- function(discriminability,
                      bias,
                      dv,
                      trial_type_var = "trial_type",
                      data,
                      distribution = "gaussian",
                      correlate_sdt_params = T,
                      # tests = "Wald",
                      control = NULL) {

  #### Check input
  if (typeof(discriminability) != "language") stop("'discriminability' must be a formula'.")
  if (typeof(bias) != "language") stop("'bias' must be a 'formula'.")

  if (typeof(dv) != "character") stop("'dv' must be of type 'character'.")
  if (is.null(data[[dv]])) stop(paste("Given dependent variable", dv, "not in data."))
  if (length(unique(data[[dv]])) != 2) stop("dv must be a binary variable.")
  if (all(sort(unique(data[[dv]])) != c(0, 1))) stop("dv must be coded as 0 ('noise' response) and 1 ('signal' response)")

  if (typeof(trial_type_var) != "character") stop("'trial_type_var' must be of type 'character'.")
  if (is.null(data[[trial_type_var]])) stop(paste("Given trialtype variable", trial_type_var, "not in data."))
  # TODO: if you have a predictor that only affects sensitivity (such as strength in the context of
  # memory, this won't work) -> maybe allow a ternary variable then (maybe with a warning)
  if (all(sort(unique(data[[trial_type_var]])) != c(-1, 1))) {
    stop("'trial_type_var' must be a numeric binary variable coding signal trials with 1 and noise trials with -1.")
  }
  if (class(data[[trial_type_var]]) != "numeric")
    stop("'trial_type_var' must be a numeric binary variable coding signal trials with 1 and noise trials with -1.")

  if (typeof(correlate_sdt_params) != "logical") stop("'correlate_sdt_params' must be of type 'logical'.")

  distribution <- standardize_dist_input(distribution)
  if (is.null(distribution)) stop("Distribution must be gaussian, logistic, or gumbel-min.")

  #### Prep & fit model
  if (distribution == "gumbel-min") {
    glmer_formula <- construct_glmer_formula(discriminability, bias, "dv_rev", mm,
                                             correlate_sdt_params = correlate_sdt_params)
  } else {
    glmer_formula <- construct_glmer_formula(discriminability, bias, dv, mm,
                                             correlate_sdt_params = correlate_sdt_params)
  }

  mm_all <- construct_modelmatrices(discriminability, bias, data, trial_type_var, distribution)
  m_frames <- mm_all[["m_frames"]]
  mm <- mm_all[["mm"]]


  # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
  # (y, ID)

  fit_obj <- fit_glmm(glmer_formula, data, mm, distribution, dv, control)
  # TODO: random effects post-processing -> for the summary method

  # Check backend stuff
  if (! is.null(summary(fit_obj)$objClass[1])) {
    print("lme4 was used to fit the model.")
    backend <- "lme4"
  } else if (any(class(fit_obj) == "glm")) {
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
