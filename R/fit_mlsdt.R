#' Convert two given formulas for the two SDT parameters to a glmer() formula
#'
#' @param formula_mu formula for sensitivity
#' @param formula_lambda formula for response bias
#' @param dv name of dependent variable
#' @param mm model matrices (corresponding to the formulas) -> only necessary for formulas with uncorrelated random effects
#' @param param_idc optional vector of parameters indices to be removed to construct a reduced formula
#' @param remove_from_mu optional argument to indicate whether the to-be-removed parameter
#'    should be removed from mu or lambda model matrix
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
construct_glmer_formula <- function(formula_mu, formula_lambda, dv, mm = NULL,param_idc = NULL, remove_from_mu = NULL) {

  # Handle missing random effects terms
  if (is.null(lme4::findbars(formula_lambda)) | is.null(lme4::findbars(formula_mu))) {
    message("At least one formula does not contain any random-effects terms. Returning NULL.")
    return()
  }

  # check if the random grouping factor is the same for mu and lambda
  rdm_fac_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][2]
  rdm_fac_lambda <- strsplit(as.character(lme4::findbars(formula_lambda)), "\\|")[[1]][2]

  if (rdm_fac_mu != rdm_fac_lambda) {
    message("Random grouping factors must be the same for sensitivity and response bias.")
    return()
  } else rdm_fac <- rdm_fac_mu


  # Handle given random effects-structure and convert to the corresponding internal formula
  # Logic: Variables are in the same parentheses if they are correlated (no "||" stuff)

  rd_pred_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })
  rd_pred_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })

  # Case 1: Everything Correlated
  # -> one random-effects term
  # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)
  # 0 to suppress automatic intercept (is contained in the mm)
  if (length(lme4::findbars(formula_lambda)) == 1 & length(lme4::findbars(formula_mu)) == 1) {
    rdm_formula <- paste("(0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | ", rdm_fac, ")", sep = "")
  } else {
    message("Given random-effects structure contains uncorrelated terms. Modeling all random effects
            parametes as uncorrelated since a mix of correlated and uncorrelated terms is not
            supported at the moment.")
    # Case 2: No Correlations
    # -> as many random-effects terms as predictors in the model
    rdm_formula_lambda <- paste(sapply(1:ncol(mm[["rdm_lambda"]]), function(x) {
      return(paste("(0 + mm[['rdm_lambda']][, ", x, "] | ", rdm_fac, ")", sep =""))
    }), collapse = "+")

    rdm_formula_mu <- paste(sapply(1:ncol(mm[["rdm_mu"]]), function(x) {
      return(paste("(0 + mm[['rdm_mu']][, ", x, "] | ", rdm_fac, ")", sep =""))
    }), collapse = "+")

    rdm_formula <- paste(rdm_formula_lambda, rdm_formula_mu, sep = " + ")

  }
    # Case 3: Mix -> TODO for later
    # -> some terms correlated, some not
    # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)

  # make the full model formula if no parameter index to be removed is given
  if (is.null(param_idc))  fixed_formula <- "0 + mm[['lambda']] + mm[['mu']]"
  else {
    # make a reduced model formula
    param_idc_char <- deparse(param_idc)
    if (grepl(":", param_idc_char)) {
      print("hi")
      param_idc_char <- paste("c(", param_idc_char, ")", sep = "")
    }
    if (remove_from_mu) {
      #print(deparse(param_idc))
      # deparse() vector for nonstandard evaluation
      term_reduced <- paste("mm[['mu']][, -", param_idc_char, ']', sep = "")
      fixed_formula <- paste("0 + mm[['lambda']] + ", term_reduced, sep = "")
    } else {
      term_reduced <- paste("mm[['lambda']][, -", param_idc_char, ']', sep = "")
      fixed_formula <- paste("0 + ", term_reduced, " + mm[['mu']]", sep = "")
    }
  }

  glmer_formula <- as.formula(paste(dv, " ~ ", fixed_formula, " + ", rdm_formula, sep = ""),
                              # parent.frame() sets scope of the parent environment (i.e., where the function
                              # is called from) for the formula -> necessary such that model matrices can be found
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
#' construct_modelmatrices(
#'   formula_mu = ~ x1 + (1 | ID)
#'   formula_lambda = ~ x1 + (1 | ID)
#'   data = data
#' )
construct_modelmatrices <- function(formula_mu,
                            formula_lambda,
                            dv,
                            data,
                            trial_type_var = "trial_type") {
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
  mm_mu <- stats::model.matrix(lme4::nobars(formula_mu), data = data)
  mm_mu <- mm_mu * trial_type_ef

  # mm_rdm_lambda
  rdm_pred_lambda <- paste(sapply(lme4::findbars(formula_lambda), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
    }), collapse = "+")

  mm_rdm_lambda <- stats::model.matrix(formula(paste("~", rdm_pred_lambda, sep = "")),
                                          data = data)

  # mm_rdm_mu
  # rdm_pred_mu <- strsplit(as.character(lme4::findbars(formula_mu)), "\\|")[[1]][1]
  rdm_pred_mu <- paste(sapply(lme4::findbars(formula_mu), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
    }), collapse = "+")

  mm_rdm_mu <- stats::model.matrix(formula(paste("~ ", rdm_pred_mu, sep = "")),
                                      data = data)
  # As above: multiply this with trial_type variable to code interaction with that factor
  mm_rdm_mu <- mm_rdm_mu * trial_type_ef


  # the model matrices consist only of the predictor variables for mu and lambda
  # for the fixed and random effects, respectively
  # via:
  # ef <- attr(terms(formula_lambda), "term.labels")
  # mapping <- attr(mm_lambda, "assign")
  # all parameters of a model term (i.e., all effect-coded predictors for a 3-level
  # factor x1) can be removed from the model matrices for LRTs of nested models

  return(list(
    "mu" = mm_mu,
    "lambda" = mm_lambda,
    "rdm_mu" = mm_rdm_mu,
    "rdm_lambda" = mm_rdm_lambda
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
  print(correlations_in_model)
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
                      backend = "lme4",
                      max_iter = 1e6) {

  if (backend == "lme4") {
    mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
    glmer_formula_full <- construct_glmer_formula(formula_mu, formula_lambda, dv, mm)

    # glmer() call consists of a mix of model matrices (model_data) and variables in "data"
    # (y, ID)

    # Do the random-effects structure selection
    # Starting with the faster optimizer (nAGQ = 0)
    glmer_formula <- glmer_formula_full
    count <- 0
    #while (TRUE) {
    #  print("Fitting")
      # 1. Fit the full model
      fit_obj_fast <- lme4::glmer(glmer_formula,
                             data = data,
                             family = binomial(link = "probit"),
                             # this is only for testing speed -> changed for actual use
                             nAGQ = 0)
      #if (isSingular(fit_obj_fast) | length(fit_obj_fast@optinfo$conv$lme4) > 0) {
    #  count <- count + 1
    #  if (count > 1) break;
    #  glmer_formula <- reduce_random_effects_structure(glmer_formula, fit_obj_fast)
    #  print(glmer_formula)

      #} else {
        # 2. Try including the correlations between random effects again

        # If the model converges, use that one, else the one before
      #  break;
      #}


    #}
    fit_obj <- fit_obj_fast
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
  } else {
    message("Only lme4 backend supported at the moment.")
    return()
  }

  # Post-Processing the lme4 output
  # TODO: this does not work for factors (rownames() is the problem)
  coefs_lambda <- summary(fit_obj)$coefficients[grepl("lambda", rownames(summary(fit_obj)$coefficients)), ]
  #rownames(coefs_lambda) <- gsub('mm', "", rownames(coefs_lambda))
  #rownames(coefs_lambda) <- colnames(mm[["lambda"]])
  if (is.null(nrow(coefs_lambda))) {
    coefs_lambda <- t(data.frame(coefs_lambda))
  } else {
    coefs_lambda <- data.frame(coefs_lambda)
  }
  rownames(coefs_lambda) <- colnames(mm[["lambda"]])

  coefs_mu <- summary(fit_obj)$coefficients[grepl("mu", rownames(summary(fit_obj)$coefficients)), ]
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


fit_submodels <- function(formula_mu, formula_lambda, dv, data, mm, type = 3, test_intercepts = F) {
  # Default behavior: generate reduced models for all _variables_ (not factors) in the formula
  # -> correspond to multiple model parameters for factors with more than two levels

  if (test_intercepts) {
    if (type == 2) {
      message("Cannot test intercepts for type II sums of squares (there's nothing left in the model).
                            Returning NULL.")
      return()
    }
    if (ncol(mm[["mu"]]) <= 1 & ncol(mm[["lambda"]]) <= 1) {
      message("Can only test intercepts if there are predictors in the model. Returning NULL.")
      return()
    }
    print(ncol(mm[["mu"]]))
    print(ncol(mm[["lambda"]]))

    range_lambda <- 0:max(attr(mm[["lambda"]], "assign"))
    range_mu <- 0:max(attr(mm[["mu"]], "assign"))
  } else {
    range_lambda <- 1:max(attr(mm[["lambda"]], "assign"))
    range_mu <- 1:max(attr(mm[["mu"]], "assign"))
  }
  print(range_lambda)
  if (type == 3) {
    reduced_formulas_lambda <- lapply(range_lambda, function(x) {
      to_remove <- which(attr(mm[["lambda"]], "assign") == x)
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = F)
      })
    reduced_formulas_mu <- lapply(range_mu, function(x) {
      to_remove <- which(attr(mm[["mu"]], "assign") == x)
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = T)
      })
  } else if (type == 2) {
    # full_formulas contain as many formulas as there are levels of interaction in the model
    # (e.g., three formulas when there are main effects, two-way and three-way interactions)
    # the last formula and model is always the full model
    # 1st element: only main effects
    # 2nd element: main effects and two-way interactions
    # etc.
    # For tests of mu predictors, the full lambda formula is included in the model (and vice versa)

    assigns <- attr(mm[["lambda"]], "assign")
    orders <- attr(terms(lme4::nobars(formula_lambda)), "order")
    # no need to fit new full models if there are only main effects in the model
    if (max(orders) < 2) full_formulas_lambda <- NULL
    else {
      n_orders_lambda <- 1:(max(orders) - 1)
      #print(n_orders_lambda)
      full_formulas_lambda <- lapply(n_orders_lambda, function(x) {
        # add 0 for the intercept (which is dropped since it has no order)
        to_remove <- which(c(0, orders[assigns]) > x)
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = F)
      })
    }
    # reduced_formulas are constructed for all predictors in the model for both mu and lambda
    # -> in order of the order of terms in the model matrix
    reduced_formulas_lambda <- lapply(range_lambda, function(x) {
      order_of_x <- orders[assigns][x]
      to_remove <- c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x))

      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = F)
      })

    assigns <- attr(mm[["mu"]], "assign")
    orders <- attr(terms(lme4::nobars(formula_mu)), "order")
    if (max(orders) < 2) full_formulas_mu <- NULL
    else {
      n_orders_mu <- 1:(max(orders) - 1)

      full_formulas_mu <- lapply(n_orders_mu, function(x) {
        to_remove <- which(c(0, orders[assigns]) > x)
        print(to_remove)
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = T)
      })
    }

    reduced_formulas_mu <- lapply(range_mu, function(x) {
      order_of_x <- orders[assigns][x]
      to_remove <- c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x))
      print(to_remove)
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = T)
      })
  }

  # set names of list to param names
  if (test_intercepts) {
    names(reduced_formulas_lambda) <- c("Intercept", paste(attr(terms(nobars(formula_lambda)), "term.labels"), "lambda", sep = "_"))
    names(reduced_formulas_mu) <- c("Intercept", paste(attr(terms(nobars(formula_mu)), "term.labels"), "mu", sep = "_"))
  } else {
    names(reduced_formulas_lambda) <- paste(attr(terms(nobars(formula_lambda)), "term.labels"), "lambda", sep = "_")
    names(reduced_formulas_mu) <- paste(attr(terms(nobars(formula_mu)), "term.labels"), "mu", sep = "_")
  }

  if (type == 2) {
    print("Fitting full models")
    full_fits_lambda <- lapply(full_formulas_lambda, function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      fit_full <- lme4::glmer(formula_tmp,
                                 data = data,
                                 family = binomial(link = "probit"),
                                 # this is only for testing speed -> changed for actual use
                                 nAGQ = 0)
    })
    full_fits_mu <- lapply(full_formulas_mu, function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      fit_full <- lme4::glmer(formula_tmp,
                              data = data,
                              family = binomial(link = "probit"),
                              # this is only for testing speed -> changed for actual use
                              nAGQ = 0)
    })
    print("Fitting reduced models")
    full_formulas_lambda_char <- sapply(full_formulas_lambda, function(x) {as.character(x)[3]})
    reduced_fits_lambda <- lapply(reduced_formulas_lambda, function(formula_tmp) {
      if (as.character(formula_tmp)[3] %in% full_formulas_lambda_char) {
        print(paste("Copying ", formula_tmp, sep = ""))
        which_model_equal <- which(full_formulas_lambda_char == as.character(formula_tmp)[3])
        fit_reduced <- full_fits_lambda[[which_model_equal]]
      } else {
        print(paste("Fitting ", formula_tmp, sep = ""))
        fit_reduced <- lme4::glmer(formula_tmp,
                                   data = data,
                                   family = binomial(link = "probit"),
                                   # this is only for testing speed -> changed for actual use
                                   nAGQ = 0)
      }
    })

    full_formulas_mu_char <- sapply(full_formulas_mu, function(x) {as.character(x)[3]})
    reduced_fits_mu <- lapply(reduced_formulas_mu, function(formula_tmp) {
      if (as.character(formula_tmp)[3] %in% full_formulas_mu_char) {
        print(paste("Copying ", formula_tmp, sep = ""))
        which_model_equal <- which(full_formulas_mu_char == as.character(formula_tmp)[3])
        print(which_model_equal)
        fit_reduced <- full_fits_mu[[which_model_equal]]
      } else {
        print(paste("Fitting ", formula_tmp, sep = ""))
        fit_reduced <- lme4::glmer(formula_tmp,
                                  data = data,
                                  family = binomial(link = "probit"),
                                  # this is only for testing speed -> changed for actual use
                                  nAGQ = 0)
      }

    })
  } else {
    # fit reduced models
    print("Fitting reduced models")
    reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      fit_reduced <- lme4::glmer(formula_tmp,
                                 data = data,
                                 family = binomial(link = "probit"),
                                 # this is only for testing speed -> changed for actual use
                                 nAGQ = 0)
    })
  }

  if (type == 2) return(list("reduced_fits_lambda" = reduced_fits_lambda,
                             "full_fits_lambda" = full_fits_lambda,
                             "reduced_fits_mu" = reduced_fits_mu,
                             "full_fits_mu" = full_fits_mu))
  else return(reduced_fits)
}


compute_LRTs <- function(fit_obj = NULL, formula_mu, formula_lambda, dv, data,
                         mm, type = 3, test_intercepts = F) {
  # only removes fixed effect, corresponding random slopes stay in the reduced model
  # So far, only type = 3 is implemented
  # TODO: also test intercept of both model matrices --> mean sensitivity and response bias

  if (! type %in% c(2, 3)) stop("Please set type to 2 or 3. Returning NULL.")

  submodels <- fit_submodels(formula_mu, formula_lambda, dv, data, mm, type, test_intercepts)
  if (type == 3) {
    reduced_fits <- submodels
    LL_full <- stats::logLik(fit_obj)
    LRT_results <- sapply(reduced_fits, function(fit_tmp) {

      LL_reduced <- stats::logLik(fit_tmp)
      chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
      p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
      df <- stats::df.residual(fit_tmp) - stats::df.residual(fit_obj)

      return(data.frame(
        # columns names inspired by afex
        "deviance_full" = -2 * LL_full,
        "deviance_reduced" = -2 * LL_reduced,
        "df.LRT" = df,
        "Chisq" = chisq,
        "p.value" = p_value
      ))
    })
  }

  else if (type == 2) {
    #return(submodels)
    reduced_fits_lambda <- submodels[["reduced_fits_lambda"]]
    reduced_fits_mu <- submodels[["reduced_fits_mu"]]
    # full fit for the highest-order effects is the already-fitted full model
    full_fits_lambda <- append(submodels[["full_fits_lambda"]], fit_obj)
    full_fits_mu <- append(submodels[["full_fits_mu"]], fit_obj)

    # assigns <- attr(mm[["lambda"]], "assign")
    orders_lambda <- attr(terms(lme4::nobars(formula_lambda)), "order")
    #print(max(orders) == length(full_fits))
    LRTs_lambda <- sapply(1:length(reduced_fits_lambda), function(fit_ind) {
      fit_tmp <- reduced_fits_lambda[[fit_ind]]
      LL_reduced <- stats::logLik(fit_tmp)
      LL_full <- stats::logLik(full_fits_lambda[[orders[fit_ind]]])
      chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
      p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
      df <- stats::df.residual(fit_tmp) - stats::df.residual(full_fits_lambda[[orders[fit_ind]]])
      return(data.frame(
        # columns names inspired by afex
        "deviance_full" = -2 * LL_full,
        "deviance_reduced" = -2 * LL_reduced,
        "df.LRT" = df,
        "Chisq" = chisq,
        "p.value" = p_value
      ))
    })


    orders_mu <- attr(terms(lme4::nobars(formula_mu)), "order")
    LRTs_mu <- sapply(1:length(reduced_fits_mu), function(fit_ind) {
      fit_tmp <- reduced_fits_mu[[fit_ind]]
      LL_reduced <- stats::logLik(fit_tmp)
      LL_full <- stats::logLik(full_fits_mu[[orders[fit_ind]]])
      chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
      p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
      df <- stats::df.residual(fit_tmp) - stats::df.residual(full_fits_mu[[orders[fit_ind]]])
      return(data.frame(
        # columns names inspired by afex
        "deviance_full" = -2 * LL_full,
        "deviance_reduced" = -2 * LL_reduced,
        "df.LRT" = df,
        "Chisq" = chisq,
        "p.value" = p_value
      ))
    })

    LRT_results <- cbind(LRTs_lambda, LRTs_mu)
  }


  # compute LRTs
  return(list(
    #"reduced_fits" = c(reduced_fits_lambda, reduced_fits_mu),
    "LRTs" = t(LRT_results)
    ))
}


