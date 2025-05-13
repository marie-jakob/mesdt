
#' Construct Model Matrices From the Given Data to give to glmer()
#'
#' @param formula_mu formula for sensitivity (mu)
#' @param formula_lambda formula for response bias (lambda)
#' @param data dataset used to construct the model data
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#'
#' @return list of model matrices (fixed and random for mu and lambda)
construct_modelmatrices <- function(formula_mu,
                                    formula_lambda,
                                    data,
                                    trial_type_var = "trial_type",
                                    distribution = "gaussian") {
  # So far: only tested for categorical predictors


  # model matrix for response bias
  # -> effects on lambda are simply main effects in the model -> predictors
  # can be included in the model without any transformation
  mm_lambda <- stats::model.matrix(lme4::nobars(formula_lambda),
                                   data = data)
  # column names of model matrix are in attr(mm_lambda, "dimnames")[[2]]

  # model matrix for sensitivity
  # -> effects on mu are interactions with trial_type variable
  # set sum contrasts for transformation of parameters later
  # -> corresponds to SDT parametrization with 0 between the two distributions
  data[["trial_type"]] <- data[[trial_type_var]]
  # contrasts(data[["trial_type"]]) <- contr.sum(2)

  # -> Intercept of model matrix becomes the mean sensitivity
  # coded with 0.5 and -0.5 such that intercept and effects can be interpreted
  # directly as increases in sensitivity
  trial_type_ef <- stats::model.matrix(~ trial_type, data = data)[, 2] * 0.5
  if (distribution == "gumbel-min") trial_type_ef <- trial_type_ef * (-1)

  mm_mu_raw <- stats::model.matrix(lme4::nobars(formula_mu), data = data)

  mm_mu <- mm_mu_raw * trial_type_ef

  to_return <- list(
    "mm" = list(
      "mu" = mm_mu,
      "lambda" = mm_lambda,
      "mu_raw" = mm_mu_raw
    )
  )

  # mm_rdm_lambda
  # check if the random grouping factor is the same for mu and lambda
  rdm_facs_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })
  rdm_facs_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })
  rdm_pred_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })
  rdm_pred_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })

  for (rdm_fac in unique(rdm_facs_mu)) {
    # get all predictors "belonging" to one random effects grouping factor
    rdm_pred_mu_tmp <- paste(rdm_pred_mu[which(rdm_facs_mu == rdm_fac)], collapse = "+")
    form_tmp <- formula(paste("~", rdm_pred_mu_tmp, sep = ""))
    mm_rdm_mu_tmp <- stats::model.matrix(form_tmp,
                                         data = data)
    mm_rdm_mu_tmp <- mm_rdm_mu_tmp * trial_type_ef
    name_rdm_mu_tmp <- paste("rdm_mu_", rdm_fac, sep = "")
    to_return[["mm"]][[name_rdm_mu_tmp]] <- mm_rdm_mu_tmp
  }
  for (rdm_fac in unique(rdm_facs_lambda)) {
    # get all predictors "belonging" to one random effects grouping factor
    rdm_pred_lambda_tmp <- paste(rdm_pred_lambda[which(rdm_facs_lambda == rdm_fac)], collapse = "+")
    #print(rdm_pred_lambda_tmp)
    form_tmp <- formula(paste("~", rdm_pred_lambda_tmp, sep = ""))
    mm_rdm_lambda_tmp <- stats::model.matrix(form_tmp,
                                             data = data)
    name_rdm_lambda_tmp <- paste("rdm_lambda_", rdm_fac, sep = "")
    to_return[["mm"]][[name_rdm_lambda_tmp]] <- mm_rdm_lambda_tmp
  }

  # Additionally create model frames to use later for emmeans
  m_frame_lambda <- stats::model.frame(lme4::nobars(formula_lambda),
                                        data = data)
  m_frame_mu <- stats::model.frame(lme4::nobars(formula_mu),
                                   data = data)
  to_return[["m_frames"]][["lambda"]] <- m_frame_lambda
  to_return[["m_frames"]][["mu"]] <- m_frame_mu

  # the model matrices consist only of the predictor variables for mu and lambda
  # for the fixed and random effects, respectively
  # via:
  # ef <- attr(terms(formula_lambda), "term.labels")
  # mapping <- attr(mm_lambda, "assign")
  # all parameters of a model term (i.e., all effect-coded predictors for a 3-level
  # factor x1) can be removed from the model matrices for LRTs of nested models

  return(to_return)
}
