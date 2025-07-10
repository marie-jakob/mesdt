#' @importFrom stats glm
#'
fit_glmm <- function(glmer_formula,
                     data,
                     mm,
                     distribution,
                     dv = NULL,
                     control = NULL) {
  # mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
  # get global options
  if (! (options("mesdt.backend") %in% c("lme4", "glmmTMB"))) {
    message("Only lme4 and glmmTMB backends are supported at the moment. Defaulting to lme4.")
    (options("mesdt.backend" = "lme4"))
  }

  if (options("mesdt.backend") == "glmmTMB" & !requireNamespace("glmmTMB", quietly = TRUE)) {
    message("Package \"glmmTMB\" must be installed and loaded to use it as backend. Setting backend to lme4.")
    options("mesdt.backend" = "lme4")
  }

  # only applies to lme4, defaults to 0 at the moment
  # -> might be changed later to a function argument
  # nAGQ <<- ifelse(is.null(options("nAGQ")$nAGQ), 0, options("nAGQ"))
  #print(nAGQ)

  # cloglog link is used for gumbel-max and gumbel-min distribution
  # -> For gumbel-min, dependent variable is reversed
  lnk_fun <- ifelse(distribution == "gaussian", "probit",
                    ifelse(distribution == "logistic", "logit", "cloglog"))

  # reverse-code dv for gumbel-min distribution
  if (distribution == "gumbel-min") {
    print("reversing dv")
    data[["dv_rev"]] <- ifelse(data[[dv]] == 0, 1, 0)
    #data[[dv]] <- ifelse(data[[dv]] == 0, 1, 0)
  }

  # Fit a GLM if there are no random effects
  #print(glmer_formula)
  if (is.null(lme4::findbars(glmer_formula))) {
    message("Formula does not contain any random terms. Fitting a single-level model.")
    fit_obj <- stats::glm(glmer_formula,
                          data = data,
                          family = binomial(link = lnk_fun))
  } else if ((options("mesdt.backend") == "lme4")) {
    if (is.null(control) |
        ! "glmerControl" %in% attr(control, "class")) {
      if (! is.null(control)) message("Invalid control argument for lme4. Using the default settings.")
      fit_obj <- lme4::glmer(glmer_formula,
                             data = data,
                             family = binomial(link = lnk_fun),
                             # this is only for testing speed -> changed for actual use
                             nAGQ = 0)
    } else {
      fit_obj <- lme4::glmer(glmer_formula,
                             data = data,
                             family = binomial(link = lnk_fun),
                             # this is only for testing speed -> changed for actual use
                             nAGQ = 0,
                             control = control)
    }

  } else if ((options("mesdt.backend") == "glmmTMB")) {
    #print("fitting with glmmTMB")
    if (is.null(control) |
        "glmerControl" %in% attr(control, "class")) {
      if (! is.null(control)) message("Invalid control argument for glmmTMB. Using the default settings.")
      fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                  data = data,
                                  family = binomial(link = lnk_fun))
    } else {
      fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                  data = data,
                                  family = binomial(link = lnk_fun),
                                  control = control)
    }
  }

  return(fit_obj)
}
