fit_glmm <- function(glmer_formula,
                      data) {
  # get global options
  backend = ifelse(options("backend") == "", "lme4", options("backend"))
  if (! backend %in% c("lme4", "glmmTMB")) {
    message("Only lme4 and glmmTMB backends are supported at the moment. Defaulting to lme4.")
    backend = "lme4"
  }
  # only applies to lme4, defaults to 0 at the moment
  # -> might be changed later to a function argument
  nAGQ = ifelse(is.null(options("nAGQ")$nAGQ), 0, options("nAGQ"))

  if (backend == "lme4") {
    fit_obj <- lme4::glmer(glmer_formula,
                           data = data,
                           family = binomial(link = "probit"),
                           # this is only for testing speed -> changed for actual use
                           nAGQ = nAGQ)
  } else if (backend == "glmmTMB") {
    fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                data = data,
                                family = binomial(link = "probit"))
  }
  return(fit_obj)
}

