#' Convert two given formulas for the two SDT parameters to an lme4 formula
#'
#' @param form_mu formula for sensitivity
#' @param form_lambda formula for response bias
#' @param dv name of dependent variable
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#' @param random name of the variable of the random factor
#' @param between predictors varying between levels of the random factor
#' @param within predictors varying within levels of the random factor
#'
#' @return lme4 formula
#' @export
#'
#' @importFrom stats formula
#'
#' @examples
#' make_glmer_formula(
#'   form_mu = mu ~ x_test,
#'   form_lambda = lambda ~ x_test,
#'   dv = "Y",
#'   trial_type_var = "trial_type",
#'   random = "ID",
#'   within = c("x_test")
#'  )
make_glmer_formula <- function(form_mu, form_lambda, dv,
                             trial_type_var = "trial_type",
                             random,
                             between = NULL,
                             within = NULL) {

  # TODO: add catches for submitted variables

  # build a new formula according to lme4 syntax
  # fixed effects
  fixed_terms <- paste(dv, "~", trial_type_var, sep = " ")

  # response bias
  if (! is.null(form_lambda)) {
    pred_lambda <- all.vars(form_lambda)[2:length(all.vars(form_lambda))]
    for (pred in pred_lambda) fixed_terms <- paste(fixed_terms, pred, sep = " + ")
  }

  # extract predictors for sensitivity and response bias from formulas
  # sensitivity
  if (! is.null(form_mu)) {
    pred_mu <- all.vars(form_mu)[2:length(all.vars(form_mu))]
    for (pred in pred_mu) {
      new_term <- paste(pred, ":", trial_type_var, sep = "")
      fixed_terms <- paste(fixed_terms, new_term, sep = " + ")
    }
  }

  # random effects
  random_terms <- paste("(", trial_type_var, sep = "")

  for (i in within) {
    if (! is.null(form_lambda)) {
      if (i %in% pred_lambda) random_terms <- paste(random_terms, i, sep = " + ")
    }
    if (! is.null(form_mu)) {
      if (i %in% pred_mu) {
        new_term <- paste(pred, ":", trial_type_var, sep = "")
        random_terms <- paste(random_terms, new_term, sep = " + ")
      }
    }
  }
  random_terms <- paste(random_terms, " | ", random, ")", sep = "")

  form_sdt <- as.formula(paste(fixed_terms, random_terms, sep = " + "),
                         env = globalenv())

  return(form_sdt)
}


#' Fit the GLMM
#'
#' @param formula The glmer model formula
#' @param data data the model should be fitted to
#' @param program whether to use R (lme4) or Julia (MixedModels)
#'
#' @return lme4 fit
#' @export
#'
#' @import lme4
#'
#' @examples
fit_glmm <- function(formula, data, program = "R") {

  if (program == "R") {
    fit <- lme4::glmer(formula, data,
                       family = binomial(link = "probit"))

  }

}



#' Title
#'
#' @param fit_obj
#' @param dat
#'
#' @return
#' @export
#'
#' @import broom.mixed
#'
#' @examples
transform_to_sdt <- function(fit_obj, dat, trial_type_var, pred_mu, pred_lambda) {
  # Start with:
  # categorical, effect-coded predictors x1 * x2 (two levels)

  # extract population-level estimates
  Beta <- broom.mixed::tidy(fit_obj, effects = "fixed")
  Beta_lambda <- Beta[! grepl("trial_type", Beta$term), ]
  Beta_mu <- Beta[grepl("trial_type", Beta$term), ]

  # make predictor design matrix
  # pred_combinations_lambda <- expand.grid()
  pred_combinations_lambda <- t(unique(dat[, pred_lambda]))
  design_matrix_lambda <- matrix(c(rep(1, ncol(pred_combinations_lambda)), pred_combinations_lambda),
                          ncol = ncol(pred_combinations_lambda), nrow = nrow(pred_combinations_lambda) + 1,
                          byrow = T)
  Lambda <- design_matrix_lambda %*% Beta_lambda$estimate * (-1)

  pred_combinations_mu <- t(unique(dat[, pred_mu]))
  design_matrix_mu <- matrix(c(rep(1, ncol(pred_combinations_mu)), pred_combinations_mu),
                                 ncol = ncol(pred_combinations_mu), nrow = nrow(pred_combinations_mu) + 1,
                                 byrow = T)
  Mu <- design_matrix_mu %*% Beta_mu$estimate * 2


  return(list(
    "Lambda" = Lambda,
    "Mu" = Mu
  ))
}
