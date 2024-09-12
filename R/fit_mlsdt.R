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
                      correlate_sdt_params = T,
                      max_iter = 10,
                      fast = T,
                      selection = NULL) {
  mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
  glmer_formula_full <- construct_glmer_formula(formula_mu, formula_lambda, dv, mm,
                                                correlate_sdt_params = correlate_sdt_params)

  # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
  # (y, ID)

  # Do the random-effects structure selection
  # Starting with the faster optimizer (nAGQ = 0)
  glmer_formula <- glmer_formula_full
  count <- 0; convergence <- F
  all_fits <- list(); all_forms <- list();
  tbr_indices <- c();
  while (TRUE) {
    if (convergence | count == max_iter) break;
    print("Hi")
    message(paste("Fitting model", as.character(glmer_formula)[2]))
    fit_tmp <- fit_glmm(glmer_formula, data, mm)

    # If no selection strategy is provided, return the fit of the full model
    if (is.null(selection)) convergence <- T
    # else reduce the model
    else if (isSingular(fit_tmp$fit_obj) | length(fit_tmp$fit_obj@optinfo$conv$lme4) > 0) {
      convergence <- F
      # 1. Remove correlations
      correlations_in_model <- ifelse(all(dim(attr(unclass(VarCorr(fit_tmp$fit_obj))$ID, "correlation")) != c(1, 1)),
                                      T, F)
      if (correlations_in_model) {
        message("Removing correlations from the model")
        glmer_formula <- construct_glmer_formula(formula_mu, formula_lambda, dv, mm,
                                                 correlate_sdt_params = correlate_sdt_params,
                                                 remove_correlations = T)
      } else {
        # Remove random terms from the model


      }

    } else convergence <- T

  }
  fit_obj <- fit_tmp
  # Continue with nAGQ = 1
  #fit_obj <- lme4::glmer(glmer_formula,
  #                       data = data,
  #                       family = binomial(link = "probit"),
  #                       # this is only for testing speed -> changed for actual use
  #                       nAGQ = 1)

  # Ensure convergence: increase iterations until the model converges
  #count = 1;
  #while (TRUE) {
  #  if (count > (max_iter / 1e4)) {
  #    stop("Final Model did not converge in the given maximal number of iterations.
  #         Please increase max_iter and try again.")
  #    break;
  #  }
  #  if (length(model_base_fin@optinfo$conv$lme4) == 0) break;

  #  print(count);
  #  pars = getME(model_base_fin, c("theta","fixef"))
  #  fit_obj <- update(fit_obj,
  #                           start = pars,
  #                           control = glmerControl(optCtrl = list(maxfun = 1e4)))
  #  print(model_base_fin@optinfo$conv$lme4)

  #  count = count + 1
  #}
  # TODO: decide whether to do the thing with adding the correlations again.

  # Post-Processing the lme4 output
  # backend = ifelse(options("backend") == "", "lme4", options("backend"))
  if (options("mlsdt.backend") == "lme4") {
    coefs_lambda <- summary(fit_obj)$coefficients[grepl("lambda", rownames(summary(fit_obj)$coefficients)), ]
    coefs_mu <- summary(fit_obj)$coefficients[grepl("mu", rownames(summary(fit_obj)$coefficients)), ]
  } else if (options("mlsdt.backend") == "glmmTMB") {
    coefs_lambda <- summary(fit_obj)$coefficients$cond[grepl("lambda", rownames(summary(fit_obj)$coefficients$cond)), ]
    coefs_mu <- summary(fit_obj)$coefficients$cond[grepl("mu", rownames(summary(fit_obj)$coefficients$cond)), ]
  }
  #rownames(coefs_lambda) <- gsub('mm', "", rownames(coefs_lambda))
  #rownames(coefs_lambda) <- colnames(mm[["lambda"]])
  if (is.null(nrow(coefs_lambda))) {
    coefs_lambda <- t(data.frame(coefs_lambda))
  } else {
    coefs_lambda <- data.frame(coefs_lambda)
  }
  rownames(coefs_lambda) <- colnames(mm[["lambda"]])

  if (is.null(nrow(coefs_mu))) {
    coefs_mu <- t(data.frame(coefs_mu))
  } else {
    coefs_mu <- data.frame(coefs_mu)
  }
  rownames(coefs_mu) <- colnames(mm[["mu"]])

  return(list(
    "fit_obj" = fit_obj,
    "Lambda" = coefs_lambda,
    "Mu" = coefs_mu
  ))
}


generate_random_terms_formula <- function(glmer_formula) {
  rdm_terms <- paste(sapply(lme4::findbars(glmer_formula), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
  }), collapse = "+")
  rdm_formula <- formula(paste("~", rdm_terms))
  return(rdm_formula)
}




reduce_random_effects_structure <- function(glmer_formula, fit_obj) {
  # If correlations are in the model, remove the correlations first
  correlations_in_model <- ifelse(all(dim(attr(unclass(VarCorr(fit_obj))$ID, "correlation")) != c(1, 1)),
                                  T, F)
  if (correlations_in_model) {
    print("Removing Correlations between random effects.")
    # remove correlations between random effects
    formula_string <- as.character(glmer_formula)
    formula_split <- strsplit(formula_string[3], "\\|")[[1]]
    new_formula_string <- paste(formula_string[2], formula_string[1], formula_split[1],
                                "||", formula_split[2])
    #print(new_formula_string)
    new_formula <- formula(new_formula_string, env = parent.frame())
  } else {
    new_formula <- NULL
  }
  # If not, check variances of the highest-order fixed effects

  # Remove the one that with the lowest variance

  return(new_formula)
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
  # ! mm_mu is different from the one in the formula function -> no
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
