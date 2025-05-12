# Try this
# @export
recover_data.mesdt_fit <- function(object, dpar = NULL, ...)  {
  if (dpar == "response bias") {
    form_tmp <- lme4::nobars(as.formula(paste(object$user_input$dv, paste(as.character(object$user_input$formula_lambda), collapse = ""))))
  } else {
    form_tmp <- lme4::nobars(as.formula(paste(object$user_input$dv, paste(as.character(object$user_input$formula_mu), collapse = ""))))
  }

  glm_tmp <- glm(form_tmp, data = object$internal$data, family = binomial("probit"))

  res <- recover_data(glm_tmp)
  return(res)
}

# Write a custom mesdt_emmeans() function, with an argument sdt_par, that
# specifies for which parameter the marginal means should be calculated

# Write an emmeans function for the mesdt_fit class, that then calls the mesdt
# function for both parameters and outputs (and saves) both
# emmeans.mesdt_fit()


# rank-deficient models abfangen!
# new argument: sdt_par = c("response_bias", "sensitivity")
# store model matrices in object
# object: whole fit object
# terms: terms of response bias formula
# xlev:
# @export
emm_basis.mesdt_fit <- function(object, trms, xlev, grid, dpar = NULL, ...) {
  if (dpar == "sensitivity") {
    m = object$internal$m_frame$mu
    contr <- attr(object$internal$mm$mu, "contrasts")
    X <- model.matrix(trms, grid, contrasts.arg = contr)
    # get coefficients for mu only

    if (object$user_input$backend == "lme4") {
      bhat <- lme4::fixef(object$fit_obj)[grep("mu", names(lme4::fixef(object$fit_obj)))]
      V <- .my.vcov(object$fit_obj)[grep("mu", rownames(stats::vcov(object$fit_obj))),
                                    grep("mu", colnames(stats::vcov(object$fit_obj)))]
    } else if (object$user_input$backend == "glmmTMB") {
      bhat <- glmmTMB::fixef(object$fit_obj)$cond[grep("mu", names(lme4::fixef(object$fit_obj)$cond))]
      V <- as.matrix(stats::vcov(object$fit_obj)$cond)[grep("mu", rownames(stats::vcov(object$fit_obj)$cond)),
                                                       grep("mu", colnames(stats::vcov(object$fit_obj)$cond))]
    }
    names(bhat) <- substr(names(bhat), 11, nchar(names(bhat)))
    rownames(V) <- substr(rownames(V), 11, nchar(rownames(V)))
    colnames(V) <- substr(colnames(V), 11, nchar(colnames(V)))
    nbasis = matrix(NA)
    # if matrix is rank deficient:
    dfargs <- list()
    dffun <- function(k, dfargs) Inf
  } else {
    m = object$internal$m_frame$lambda
    contr <- attr(object$internal$mm$lambda, "contrasts")
    X <- model.matrix(trms, grid, contrasts.arg = contr)
    #X = object$internal$mm$lambda
    # get coefficients for lambda only
    if (object$user_input$backend == "lme4") {
      bhat <- lme4::fixef(object$fit_obj)[grep("lambda", names(lme4::fixef(object$fit_obj)))]
      V <- .my.vcov(object$fit_obj)[grep("lambda", rownames(stats::vcov(object$fit_obj))),
                                    grep("lambda", colnames(stats::vcov(object$fit_obj)))]
    } else if (object$user_input$backend == "glmmTMB") {
      bhat <- glmmTMB::fixef(object$fit_obj)$cond[grep("lambda", names(glmmTMB::fixef(object$fit_obj)$cond))]
      V <- as.matrix(stats::vcov(object$fit_obj)$cond)[grep("lambda", rownames(stats::vcov(object$fit_obj)$cond)),
                                    grep("lambda", colnames(stats::vcov(object$fit_obj)$cond))]
    }
    names(bhat) <- substr(names(bhat), 15, nchar(names(bhat)))
    # only select rows and columns relevant for sensitivity

    rownames(V) <- substr(rownames(V), 15, nchar(rownames(V)))
    colnames(V) <- substr(colnames(V), 15, nchar(colnames(V)))
    nbasis = matrix(NA)
    # if matrix is rank deficient:
    dfargs <- list()
    dffun <- function(k, dfargs) Inf
  }
  return(list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs))
}

