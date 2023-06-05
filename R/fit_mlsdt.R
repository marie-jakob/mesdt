#' Convert two given formulas for the two SDT parameters to a glmer() formula
#'
#' @param formula_mu formula for sensitivity
#' @param formula_lambda formula for response bias
#' @param dv name of dependent variable
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#' @param param_idc optional vector of parameters indices to be removed to construct a reduced formula
#' @param remove_from_mu optional argument to indicate whether the to-be-removed parameter
#'    should be removed from mu or lambda modeldata
#'
#' @return lme4 formula
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
construct_glmer_formula <- function(formula_mu, formula_lambda, dv, param_idc = NULL, remove_from_mu = NULL) {

  # check if the random grouping factor is the same for mu and lambda
  random_fac_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][2]
  random_fac_lambda <- strsplit(as.character(lme4::findbars(formula_lambda)), "\\|")[[1]][2]

  if (random_fac_mu != random_fac_lambda) {
    message("Random grouping factors must be the same for sensitivity and response bias.")
    return()
  } else random_fac <- random_fac_mu


  # Handle
  if (is.null(lme4::findbars(formula_lambda)) | is.null(lme4::findbars(formula_mu))) {
    message("At least one formula does not contain any random-effects terms. Returning NULL")
    return()
  }

  # Handle given random effects-structure and convert to the corresponding internal formula
  # Logic: Variables are in the same parentheses iff they are correlated (no "||" stuff)

  # Case 1: No Correlations
  # -> as many random-effects terms as predictors in the model
  rd_pred_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
    })
  rd_pred_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
    })

  if (length(rd_pred_lambda) == length(lme4::findbars(formula_lambda)) &
      length(rd_pred_mu) == length(lme4::findbars(formula_mu))) {
    random_formula <- paste("(0 + modeldata[['random_lambda']] + modeldata[['random_mu']] || ", random_fac, ")", sep = "")
  }


  # Case 2: Everything Correlated
  # -> one random-effects term
  # -> ~ ... + (modeldata[["random_mu"]] + modeldata[["random_lambda]] | ID)
  if (length(lme4::findbars(formula_lambda)) == 1 & length(lme4::findbars(formula_mu)) == 1) {
    # 0 to suppress automatic intercept (is contained in the modeldata)
    random_formula <- paste("(0 + modeldata[['random_lambda']] + modeldata[['random_mu']] | ", random_fac, ")", sep = "")
  }


  # Case 3: Mix


  # make the full model formula if no parameter index to be removed is given
  if (is.null(param_idc))  fixed_formula <- "0 + modeldata[['lambda']] + modeldata[['mu']]"
  else {
    # make a reduced model formula
    if (remove_from_mu) {
      # deparse() vector for nonstandard evaluation
      term_reduced <- paste("modeldata[['mu']][, -", deparse(param_idc), ']', sep = "")
      fixed_formula <- paste("0 + modeldata[['lambda']] + ", term_reduced, sep = "")
    } else {
      term_reduced <- paste("modeldata[['lambda']][, -", deparse(param_idc), ']', sep = "")
      fixed_formula <- paste("0 + ", term_reduced, " + modeldata[['mu']]", sep = "")
    }
  }

  glmer_formula <- as.formula(paste(dv, " ~ ", fixed_formula, " + ", random_formula, sep = ""),
                              # parent.frame() sets scope of the parent environment (i.e., where the function
                              # is called from) for the formula -> necessary such that modeldata can be found
                              env = parent.frame())
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
  # coded with 0.5 and -0.5 such that intercept and effects can be interpreted
  # directly as increases in sensitivity
  trial_type_ef <- stats::model.matrix(~ trial_type, data = data)[, 2] * 0.5
  modeldata_mu <- stats::model.matrix(lme4::nobars(formula_mu), data = data)
  modeldata_mu <- modeldata_mu * trial_type_ef

  # modeldata_random_lambda
  random_pred_lambda <- paste(sapply(lme4::findbars(formula_lambda), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
    }), collapse = "+")

  modeldata_random_lambda <- stats::model.matrix(formula(paste("~", random_pred_lambda, sep = "")),
                                          data = data)

  # modeldata_random_mu
  # random_pred_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][1]
  random_pred_mu <- paste(sapply(lme4::findbars(formula_mu), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
    }), collapse = "+")

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
    "mu" = modeldata_mu,
    "lambda" = modeldata_lambda,
    "random_mu" = modeldata_random_mu,
    "random_lambda" = modeldata_random_lambda
  ))
}


#' Fit a multilevel signal detection theory model
#'
#' @param formula_mu Formula specifying fixed and random effects on sensitivity
#' @param formula_lambda Formula specifying fixed and random effects on response bias
#' @param dv name of the (binary) dependent variable
#' @param trial_type_var name of the variable coding signal vs. noise trials
#' @param data dataset
#' @param backend package / library used to fit the model
#'    -> only supports R lme4 (backend = "lme4" at the moment but support for Julia
#'    MixedModels and glmmTMB is planned)
#'
#' @return TODO
#' @import lme4
#' @export
#'
#' @examples
fit_mlsdt <- function(formula_mu,
                      formula_lambda,
                      dv,
                      trial_type_var = "trial_type",
                      data,
                      backend = "lme4") {

  if (backend == "lme4") {
    glmer_formula <- construct_glmer_formula(formula_mu, formula_lambda, dv)
    modeldata <- construct_modeldata(formula_mu, formula_lambda, dv, data, trial_type_var)

    # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
    # (y, ID)
    fit_obj <- lme4::glmer(glmer_formula,
                           data = data,
                           family = binomial(link = "probit"),
                           # this is only for testing speed -> changed for actual use
                           nAGQ = 0)
  } else {
    message("Only lme4 backend supported at the moment.")
    return()
  }

  coefs_lambda <- summary(fit_obj)$coefficients[grepl("lambda", rownames(summary(fit_obj)$coefficients)), ]
  rownames(coefs_lambda) <- gsub('modeldata', "", rownames(coefs_lambda))
  coefs_mu <- summary(fit_obj)$coefficients[grepl("mu", rownames(summary(fit_obj)$coefficients)), ]
  rownames(coefs_mu) <- gsub('modeldata', "", rownames(coefs_mu))

  return(list(
    "fit_obj" = fit_obj,
    "Lambda" = coefs_lambda,
    "Mu" = coefs_mu
  ))
}




#' Title
#'
#' @param fit_obj
#' @param data
#'
#' @return
#'
#' TODO: test this...
#' @examples
transform_to_sdt <- function(fit_obj, formula_lambda, formula_mu, data, level = 0.95) {

  # Get Betas
  fixef_lambda <- fixef(fit_obj)[grepl("lambda", names(fixef(fit_obj)))]
  fixef_mu <- fixef(fit_obj)[grepl("mu", names(fixef(fit_obj)))]

  # Transform response bias parameters
  coding_lambda <- unique(stats::model.matrix(lme4::nobars(formula_lambda),
                                              data = data))

  order_lambda <- attr(terms(lme4::nobars(formula_lambda)), "order")
  mapping_lambda <- attr(stats::model.matrix(lme4::nobars(formula_lambda),
                                             data = data), "assign")

  # + 1 because of the intercept
  factors_lambda <- coding_lambda[, which(order_lambda[mapping_lambda] == 1) + 1]
  est_lambda <- coding_lambda %*% fixef_lambda * (-1)

  est_lambda <- data.frame(cbind(factors_lambda, est_lambda))
  names(est_lambda) <- c(names(est_lambda)[1:ncol(est_lambda) - 1], "estimate")

  # Transform sensitivity parameters
  # ! modeldata_mu is different from the one in the formula function -> no
  # interaction with "trial_type" variable
  coding_mu <- unique(stats::model.matrix(lme4::nobars(formula_mu),
                                   data = data), data = test_data)

  order_mu <- attr(terms(lme4::nobars(formula_mu)), "order")
  mapping_mu <- attr(stats::model.matrix(lme4::nobars(formula_mu),
                                  data = data), "assign")

  # + 1 because of the intercept
  factors_mu <- coding_mu[, which(order_mu[mapping_mu] == 1) + 1]
  est_mu <- coding_mu %*% fixef_mu * 2

  est_mu <- data.frame(cbind(factors_mu, est_mu))
  names(est_mu) <- c(names(est_mu)[1:ncol(est_mu) - 1], "estimate")

  #----------------------------------------------------------------------------#
  #### Compute Wald SEs and CIs ####

  which_lambda <- which(grepl("lambda", names(fixef(fit_obj))))

  inv_hess_lambda <- as.matrix(vcov(fit_obj))[which_lambda, which_lambda]

  # SEs for parameter combinations (i.e., SDT parameters)
  # -> sqrt() of sum of elements of the inv-Hessian corresponding to added Beta weights
  # (i.e., non-zero cells of the parameter vector)

  SEs_lambda <- apply(coding_lambda, 1, function(x) {
    sqrt(sum(inv_hess_lambda[which(x != 0), which(x != 0)]))
  })

  est_lambda$SE <- SEs_lambda

  which_mu <- which(grepl("mu", names(fixef(fit_obj))))
  inv_hess_mu <- as.matrix(vcov(fit_obj))[which_mu, which_mu]

  SEs_mu <- apply(coding_mu, 1, function(x) {
    sqrt(sum(inv_hess_mu[which(x != 0), which(x != 0)]))
  })

  est_mu$SE <- SEs_mu

  # CIs
  print(est_lambda$estimate + qnorm((1 - level) / 2))
  est_lambda$CI_lower <- est_lambda$estimate + qnorm((1 - level) / 2) * est_lambda$SE
  est_lambda$CI_upper <- est_lambda$estimate + qnorm(level + ((1 - level) / 2)) * est_lambda$SE

  est_mu$CI_lower <- est_mu$estimate + qnorm((1 - level) / 2) * est_mu$SE
  est_mu$CI_upper <- est_mu$estimate + qnorm(level + ((1 - level) / 2)) * est_mu$SE

  return(list(
    "lambda_estimates" = est_lambda,
    "mu_estimates" = est_mu
  ))

}


fit_submodels <- function(formula_mu, formula_lambda, dv, data, modeldata, type = 3, test_intercepts = F) {
  # Default behavior: generate reduced models for all _variables_ (not factors) in the formula
  # -> correspond to multiple model parameters for factors with more than two levels
  reduced_formulas_lambda <- lapply(1:(ncol(modeldata[["lambda"]]) - 1), function(x) {
    to_remove <- which(attr(modeldata[["lambda"]], "assign") == x)
    construct_glmer_formula(formula_mu, formula_lambda, dv = dv,
                            param_idc = to_remove, remove_from_mu = F)
    }
  )
  #if (! test_intercepts) reduced_formulas_lambda <- reduced_formulas_lambda

  # set names of list to param names
  names(reduced_formulas_lambda) <- paste(attr(terms(nobars(formula_lambda)), "term.labels"), "lambda", sep = "_")

  reduced_formulas_mu <- lapply(1:(ncol(modeldata[["mu"]]) - 1), function(x) {
    to_remove <- which(attr(modeldata[["mu"]], "assign") == x)
    construct_glmer_formula(formula_mu, formula_lambda, dv = dv,
                            param_idc = to_remove, remove_from_mu = T)
    }
  )
  names(reduced_formulas_mu) <- paste(attr(terms(nobars(formula_mu)), "term.labels"), "mu", sep = "_")

  # fit reduced models
  reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), function(formula_tmp) {
    print(paste("Fitting ", formula_tmp, sep = ""))
    fit_reduced <- lme4::glmer(formula_tmp,
                               data = data,
                               family = binomial(link = "probit"),
                               # this is only for testing speed -> changed for actual use
                               nAGQ = 0)
    }
  )
  return(reduced_fits)
}


compute_LRTs <- function(fit_obj = NULL, formula_mu, formula_lambda, dv, data,
                         modeldata, type = 3, test_intercepts = F) {
  # only removes fixed effect, corresponding random slopes stay in the reduced model
  # So far, only type = 3 is implemented
  # TODO: how does type 2 make sense here (sensitivity params as interactions)
  # TODO: also test intercept of both model matrices --> mean sensitivity and response bias

  reduced_fits <- fit_submodels(formula_mu, formula_lambda, dv, data, modeldata, type, test_intercepts)

  LL_full <- stats::logLik(fit_obj)
  LL_reduced <- sapply(reduced_fits, stats::logLik)

  LRT_results <- sapply(reduced_fits, function(fit_tmp) {
    LL_reduced <- stats::logLik(fit_tmp)
    chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
    p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
    df <- stats::df.residual(fit_tmp) - stats::df.residual(fit_obj)

    return(data.frame(
      # columns names inspired by afex
      "deviance_full" = -2 * LL_full,
      "deviance_reduced" = -2 * LL_reduced,
      # df changes for type II-like LRTs
      "df.LRT" = df,
      "Chisq" = chisq,
      "p.value" = p_value
    ))
  })


  # in the future --> parallelize this

  # compute LRTs
  return(list(
    "reduced_fits" = reduced_fits,
    "LRTs" = t(LRT_results)
    ))
}


