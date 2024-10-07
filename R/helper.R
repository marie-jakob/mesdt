.onLoad <- function(libname, pkgname) {
  options(mlsdt.backend = "lme4")
}
fit_glmm <- function(glmer_formula,
                     data,
                     mm) {
  # mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
  # get global options

  if (! (options("mlsdt.backend") %in% c("lme4", "glmmTMB"))) {
    message("Only lme4 and glmmTMB backends are supported at the moment. Defaulting to lme4.")
    (options("mlsdt.backend" = "lme4"))
  }

  if (options("mlsdt.backend") == "glmmTMB" & !requireNamespace("glmmTMB", quietly = TRUE)) {
    message("Package \"glmmTMB\" must be installed to use it as backend. Setting backend to lme4.")
    options("mlsdt.backend" = "lme4")
  }

  # only applies to lme4, defaults to 0 at the moment
  # -> might be changed later to a function argument
  nAGQ = ifelse(is.null(options("nAGQ")$nAGQ), 0, options("nAGQ"))
  #print(nAGQ)


  # Fit a GLM if there are no random effects
  if (is.null(lme4::findbars(glmer_formula))) {
    message("Formula does not contain any random terms. Fitting a single-level model.")
    fit_obj <- stats::glm(glmer_formula,
                          data = data,
                          family = binomial(link = "probit"))
  } else if ((options("mlsdt.backend") == "lme4")) {
    fit_obj <- lme4::glmer(glmer_formula,
                           data = data,
                           family = binomial(link = "probit"),
                           # this is only for testing speed -> changed for actual use
                           nAGQ = nAGQ)
  } else if ((options("mlsdt.backend") == "glmmTMB")) {
    #print("fitting with glmmTMB")
    # mm <- construct_modelmatrices(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = data)
    fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                data = data,
                                family = binomial(link = "probit"))
  }
  return(fit_obj)
}


