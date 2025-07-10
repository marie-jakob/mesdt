#' @importFrom emmeans emm_basis recover_data
#' @method recover_data mesdt_fit
#' @export
recover_data.mesdt_fit <- function(object, dpar = NULL, ...)  {
  if (dpar == "response bias") {
    form_tmp <- lme4::nobars(as.formula(paste(object$user_input$dv, paste(as.character(object$user_input$bias), collapse = ""))))
  } else {
    form_tmp <- lme4::nobars(as.formula(paste(object$user_input$dv, paste(as.character(object$user_input$discriminability), collapse = ""))))
  }

  glm_tmp <- glm(form_tmp, data = object$internal$data, family = binomial("probit"))

  res <- recover_data(glm_tmp)
  return(res)
}


# rank-deficient models abfangen!
# new argument: sdt_par = c("response_bias", "sensitivity")
# store model matrices in object
# object: whole fit object
# terms: terms of response bias formula
#' @importFrom emmeans emm_basis recover_data
#' @importFrom emmeans .my.vcov
#' @importFrom stats model.matrix
#' @method emm_basis mesdt_fit
#' @export
emm_basis.mesdt_fit <- function(object, trms, xlev, grid, dpar = NULL, ...) {
  dpar <- match.arg(dpar, c("response bias", "discriminability", "sensitivity"))
  sdt_par = ifelse(dpar == "discriminability" |
                     dpar == "sensitivity", "mu", "lambda")
  m = object$internal$m_frame[[sdt_par]]
  contr <- attr(object$internal$mm[[sdt_par]], "contrasts")
  if (object$user_input$distribution != "gumbel-min") {
    mult <- ifelse(sdt_par == "lambda", -1, 1)
  } else {
    mult <- 1
  }
  X <- stats::model.matrix(trms, grid, contrasts.arg = contr)
  if (object$internal$backend == "lme4") {
    bhat <- lme4::fixef(object$fit_obj)[grep(sdt_par, names(lme4::fixef(object$fit_obj)))] * mult
      V <- as.matrix(.my.vcov(object$fit_obj)[grep(sdt_par, rownames(stats::vcov(object$fit_obj))),
                                    grep(sdt_par, colnames(stats::vcov(object$fit_obj)))])
  } else if (object$internal$backend == "glmmTMB") {
    bhat <- glmmTMB::fixef(object$fit_obj)$cond[grep(sdt_par, names(lme4::fixef(object$fit_obj)$cond))] * mult
      V <- as.matrix(stats::vcov(object$fit_obj)$cond)[grep(sdt_par, rownames(stats::vcov(object$fit_obj)$cond)),
                                                       grep(sdt_par, colnames(stats::vcov(object$fit_obj)$cond))]
  } else if (object$internal$backend == "glm") {
    bhat <- stats::coef(object$fit_obj)[grep(sdt_par, names(stats::coef(object$fit_obj)))] * mult
    V <- as.matrix(emmeans::.my.vcov(object$fit_obj)[grep(sdt_par, rownames(stats::vcov(object$fit_obj))),
                                  grep(sdt_par, colnames(stats::vcov(object$fit_obj)))])
  }
  sdt_n_char <- ifelse(dpar == "sensitivity", 11, 15)
  names(bhat) <- substr(names(bhat), sdt_n_char, nchar(names(bhat)))
  rownames(V) <- substr(rownames(V), sdt_n_char, nchar(rownames(V)))
  colnames(V) <- substr(colnames(V), sdt_n_char, nchar(colnames(V)))
  nbasis = matrix(NA)
  # if matrix is rank deficient:
  dfargs <- list()
  dffun <- function(k, dfargs) Inf
  return(list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs))
}

