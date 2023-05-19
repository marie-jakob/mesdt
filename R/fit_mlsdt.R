#' Convert two given formulas for the two SDT parameters to a glmer() formula
#'
#' @param formula_mu formula for sensitivity
#' @param formula_lambda formula for response bias
#' @param dv name of dependent variable
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#'
#' @return lme4 formula
#' @export
#'
#' @importFrom stats formula
#' @importFrom lme4 findbars
#'
#' @examples
#' construct_glmer_formula(
#'   formula_mu = ~ x1 + (1 | ID)
#'   formula_lambda = ~ x2 + (1 | ID)
#'   dv = "y",
#'   trial_type_var = "trial_type"
#' )
construct_glmer_formula <- function(formula_mu, formula_lambda, dv) {

  # check if the random grouping factor is the same for mu and lambda
  random_fac_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][2]
  random_fac_lambda <- strsplit(as.character(lme4::findbars(formula_lambda)), "\\|")[[1]][2]

  if (random_fac_mu != random_fac_lambda) {
    message("Random grouping factors must be the same for sensitivity and response bias.")
    return()
  } else random_fac <- random_fac_mu

  # 0 to suppress automatic intercept (is contained in the modeldata)
  random_formula <- paste("(0 + modeldata_random_lambda + modeldata_random_mu | ", random_fac, ")", sep = "")

  glmer_formula <- formula(paste(dv, "~ 0 + modeldata_lambda + modeldata_mu + ", random_formula, sep = ""))

  return(glmer_formula)
}



#' Construct Model Matrices From the Given Data to give to glmer()
#'
#' @param formula_mu formula for sensitivity (mu)
#' @param formula_lambda formula for response bias (lambda)
#' @param dv name of dependent variable
#' @param data dataset used to construct the model data
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#'
#' @return list of 4 model matrices (fixed and random for mu and lambda)
#' @export
#'
#' @importFrom lme4 nobars
#' @importFrom lme4 findbars
#' @importFrom stats formula
#' @importFrom stats model.matrix
#'
#' @examples
#' construct_modeldata(
#'   formula_mu = ~ x1 + (1 | ID)
#'   formula_lambda = ~ x1 + (1 | ID)
#'   data = data
#' )
construct_modeldata <- function(formula_mu,
                            formula_lambda,
                            dv,
                            data,
                            trial_type_var = "trial_type") {
  # So far: only tested for categorical predictors


  # modeldata for response bias
  # -> effects on lambda are simply main effects in the model -> predictors
  # can be included in the model without any transformation
  modeldata_lambda <- stats::model.matrix(lme4::nobars(formula_lambda),
                                          data = data)
  # column names of modeldata matrix are in attr(modeldata_lambda, "dimnames")[[2]]

  # modeldata for sensitivity
  # -> effects on mu are interactions with trial_type variable
  # set sum contrasts for transformation of parameters later
  # -> corresponds to SDT parametrization with 0 between the two distributions
  data[["trial_type"]] <- data[[trial_type_var]]
  contrasts(data[["trial_type"]]) <- contr.sum(2)

  # -> Intercept of modeldata matrix becomes the mean sensitivity
  trial_type_ef <- stats::model.matrix(~ trial_type, data = data)[, 2]
  modeldata_mu <- stats::model.matrix(lme4::nobars(formula_mu), data = data)
  modeldata_mu <- modeldata_mu * trial_type_ef

  # modeldata_random_lambda
  random_pred_lambda <- strsplit(as.character(lme4::findbars(formula_lambda)), "\\|")[[1]][1]

  modeldata_random_lambda <- stats::model.matrix(formula(paste("~", random_pred_lambda, sep = "")),
                                          data = data)

  # modeldata_random_mu
  random_pred_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][1]

  modeldata_random_mu <- stats::model.matrix(formula(paste("~ ", random_pred_mu, sep = "")),
                                      data = data)
  # As above: multiply this with trial_type variable to code interaction with that factor
  modeldata_random_mu <- modeldata_random_mu * trial_type_ef


  # the modeldata matrices consist only of the predictor variables for mu and lambda
  # for the fixed and random effects, respectively
  # via:
  # ef <- attr(terms(formula_lambda), "term.labels")
  # mapping <- attr(modeldata_lambda, "assign")
  # all parameters of a model term (i.e., all effect-coded predictors for a 3-level
  # factor x1) can be removed from the modeldata matrices for LRTs of nested models

  return(list(
    "modeldata_mu" = modeldata_mu,
    "modeldata_lambda" = modeldata_lambda,
    "modeldata_random_mu" = modeldata_random_mu,
    "modeldata_random_lambda" = modeldata_random_lambda
  ))
}



#' Fit the GLMM
#'
#' @param formula The glmer model formula
#' @param data data the model should be fitted to
#' @param backend whether to use R (lme4) or Julia (MixedModels)
#'
#' @return lme4 fit
#' @export
#'
#' @import lme4
#'
#' @examples
fit_glmm <- function(formula, data, backend = "lme4") {

  if (program == "R") {
    fit <- lme4::glmer(formula, data,
                       family = binomial(link = "probit"))

  }
  return()

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
