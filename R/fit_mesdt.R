#' Fit a multilevel signal detection theory model
#'
#' @param formula_mu Formula specifying fixed and random effects on sensitivity
#' @param formula_lambda Formula specifying fixed and random effects on response bias
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
fit_mesdt <- function(formula_mu,
                      formula_lambda,
                      dv,
                      trial_type_var = "trial_type",
                      data,
                      correlate_sdt_params = T,
                      tests = "Wald",
                      control = NULL) {

  if (! tests %in% c("Wald", "LRT", "Boot", "wald", "lrt", "boot")) {
    stop(paste("Tests of type ", tests, " not available. Please use Wald, LRT or boot.", sep = ""))
  }

  mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
  glmer_formula <- construct_glmer_formula(formula_mu, formula_lambda, dv, mm,
                                                correlate_sdt_params = correlate_sdt_params)

  # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
  # (y, ID)

  fit_obj <- fit_glmm(glmer_formula, data, mm, control)
  # TODO: random effects post-processing -> for the summary method

  # Compute tests
  # TODO: do I want this

  #if (tests == "Wald") {
  #  # Post-Processing the lme4 output
  #  # backend = ifelse(options("backend") == "", "lme4", options("backend"))
  #  if (options("mlsdt.backend") == "lme4") {
  #    coefs_lambda <- summary(fit_obj)$coefficients[grepl("lambda", rownames(summary(fit_obj)$coefficients)), ]
  #    coefs_mu <- summary(fit_obj)$coefficients[grepl("mu", rownames(summary(fit_obj)$coefficients)), ]
  #  } else if (options("mlsdt.backend") == "glmmTMB") {
  #    coefs_lambda <- summary(fit_obj)$coefficients$cond[grepl("lambda", rownames(summary(fit_obj)$coefficients$cond)), ]
  #    coefs_mu <- summary(fit_obj)$coefficients$cond[grepl("mu", rownames(summary(fit_obj)$coefficients$cond)), ]
  #  }
  #  #rownames(coefs_lambda) <- gsub('mm', "", rownames(coefs_lambda))
  #  #rownames(coefs_lambda) <- colnames(mm[["lambda"]])
  #  if (is.null(nrow(coefs_lambda))) {
  #    coefs_lambda <- t(data.frame(coefs_lambda))
  #  } else {
  #    coefs_lambda <- data.frame(coefs_lambda)
  #  }
  #  rownames(coefs_lambda) <- colnames(mm[["lambda"]])
#
  #  if (is.null(nrow(coefs_mu))) {
  #    coefs_mu <- t(data.frame(coefs_mu))
  #  } else {
  #    coefs_mu <- data.frame(coefs_mu)
  #  }
  #  rownames(coefs_mu) <- colnames(mm[["mu"]])
  #}

  # Check backend stuff
  if (! is.null(summary(fit_obj)$objClass[1])) {
    print("lme4 was used to fit the model.")
    backend <- "lme4"
  } else {
    print("glmmTMB was used to fit the model.")
    backend <- "glmmTMB"
  }


  #obj <- new_mesdt_fit(list(
  #  "fit_obj" = fit_obj,
    #"Lambda" = coefs_lambda,
    #"Mu" = coefs_mu,
  #  "formula_mu" = formula_mu,
  #  "formula_lambda" = formula_lambda,
  #  "dv" = dv,
  #  "trial_type_var" = trial_type_var,
  #  "backend" = backend,
  #  "correlate_sdt_params" = correlate_sdt_params
  #))

  obj <- list(
    "fit_obj" = fit_obj,
    "formula_mu" = formula_mu,
    "formula_lambda" = formula_lambda,
    "dv" = dv,
    "trial_type_var" = trial_type_var,
    "backend" = backend,
    "correlate_sdt_params" = correlate_sdt_params
    )

  return(obj)
}
