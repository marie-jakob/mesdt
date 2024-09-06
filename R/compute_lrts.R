fit_submodels <- function(formula_mu, formula_lambda, dv, data, mm, type = 3, test_intercepts = F,
                          rem_ran_ef = F) {
  # Default behavior: generate reduced models for all _variables_ (not factors) in the formula
  # -> correspond to multiple model parameters for factors with more than two levels

  if (test_intercepts) {
    #if (type == 2) {
    #  message("Cannot test intercepts for type II sums of squares (there's nothing left in the model).
    #                        Returning NULL.")
    #  return()
    #}
    #if (ncol(mm[["mu"]]) <= 1 & ncol(mm[["lambda"]]) <= 1) {
    #  message("Can only test intercepts if there are predictors in the model. Returning NULL.")
    #  return()
    #}
    if (rem_ran_ef) {
      #range_lambda and range_mu are now list (to accommodate multiple random effects)
      range_lambda <- list(); range_mu <- list();
      for (i in 1:length(mm)) {
        if (grepl("rdm_lambda", names(mm)[i])) {
          range_lambda[[names(mm)[i]]] <- 0:max(attr(mm[[names(mm)[i]]], "assign"))
        } else if (grepl("rdm_mu", names(mm)[i])) {
          range_mu[[names(mm)[i]]] <- 0:max(attr(mm[[names(mm)[i]]], "assign"))
        }
      }
    } else {
      range_lambda <- list(0:max(attr(mm[["lambda"]], "assign"))); names(range_lambda) <- c("lambda")
      range_mu <- list(0:max(attr(mm[["mu"]], "assign"))); names(range_mu) <- c("mu")
      }
  } else {
    # If there are no predictors in the model and ! test_intercepts, throw a message
    if (length(attr(mm[["lambda"]], "assign")) == 1) range_lambda <- c()
    else  range_lambda <- list(1:max(attr(mm[["lambda"]], "assign")))
    if (length(attr(mm[["mu"]], "assign")) == 1) range_mu <- c()
    else range_mu <- list(1:max(attr(mm[["mu"]], "assign")))
    if (length(range_lambda) > 0) names(range_lambda) <- "lambda";
    if (length(range_mu) > 0) names(range_mu) <- "mu";
    if (length(range_lambda) == 0 & length(range_mu) == 0) {
      message("Nothing to test in the model. Returning NULL.")
      return(NULL)
    }
  }

  if (type == 3) {
    reduced_formulas_lambda <- list()
    for (i in 1:length(range_lambda)) {
      forms_tmp <- lapply(range_lambda[[i]], function(x) {
        to_remove <- as.numeric(which(attr(mm[[names(range_lambda)[i]]], "assign") == x))
        if (length(to_remove) == dim(mm[[names(range_lambda)[i]]])[2]) to_remove = Inf
        remove_from_rdm_tmp <- ifelse(grepl("rdm", names(range_lambda)[i]), gsub("rdm_lambda_", "", names(range_lambda)[i]), "")
        return(construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                       param_idc = to_remove, remove_from_mu = F,
                                       remove_from_rdm = remove_from_rdm_tmp))
      })
      reduced_formulas_lambda <- append(reduced_formulas_lambda, forms_tmp)
    }
    reduced_formulas_mu <- list()
    for (i in 1:length(range_mu)) {
      forms_tmp <- lapply(range_mu[[i]], function(x) {
        to_remove <- as.numeric(which(attr(mm[[names(range_mu)[i]]], "assign") == x))
        if (length(to_remove) == dim(mm[[names(range_mu)[i]]])[2]) to_remove = Inf
        remove_from_rdm_tmp <- ifelse(grepl("rdm", names(range_mu)[i]), gsub("rdm_mu_", "", names(range_mu)[i]), "")
        return(construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                       param_idc = to_remove, remove_from_mu = T,
                                       remove_from_rdm = remove_from_rdm_tmp))

      })
      reduced_formulas_mu <- append(reduced_formulas_mu, forms_tmp)
    }
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
    #if (test_intercepts) orders <- c(0, orders)
    # no need to fit new full models if there are only main effects in the model
    # (if ! test_intercepts) or if there are no predictors in the model (if test_intercepts)
    # if ((max(orders) < 2) & ! test_intercepts) full_formulas_lambda <- NULL
    if ((max(c(0, orders)) + as.numeric(test_intercepts)) < 2) full_formulas_lambda <- NULL
    else {
      if (test_intercepts) n_orders_lambda <- 0:(max(orders) - 1)
      else n_orders_lambda <- 1:(max(orders) - 1)
      full_formulas_lambda <- lapply(n_orders_lambda, function(x) {
        # add 0 for the intercept (which is dropped since it has no order)
        to_remove <- as.numeric(which(c(0, orders[assigns]) > x))
        if (length(to_remove) == dim(mm[["lambda"]])[2]) to_remove = Inf
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = F)
      })
    }
    # reduced_formulas are constructed for all predictors in the model for both mu and lambda
    # -> in order of the order of terms in the model matrix
    reduced_formulas_lambda <- lapply(range_lambda[[1]], function(x) {
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- as.numeric(c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x)))
      if (length(to_remove) == dim(mm[["lambda"]])[2]) to_remove = Inf
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = F)
    })

    assigns <- attr(mm[["mu"]], "assign")
    orders <- attr(terms(lme4::nobars(formula_mu)), "order")
    #if ((max(orders) < 2) & ! test_intercepts) full_formulas_mu <- NULL
    if ((max(c(0, orders)) + as.numeric(test_intercepts)) < 2) full_formulas_mu <- NULL
    else {
      if (test_intercepts) n_orders_mu <- 0:(max(orders) - 1)
      else n_orders_mu <- 1:(max(orders) - 1)

      full_formulas_mu <- lapply(n_orders_mu, function(x) {
        to_remove <- as.numeric(which(c(0, orders[assigns]) > x))
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = T)
      })
    }

    reduced_formulas_mu <- lapply(range_mu[[1]], function(x) {
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- as.numeric(c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x)))
      if (length(to_remove) == dim(mm[["mu"]])[2]) to_remove = Inf
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = T)
    })
  }

  # set names of list to param names
  names_lambda <- c(); names_mu <- c();
  # names(reduced_formulas_lambda) <- c(); names(reduced_formulas_mu) <- c()
  if (test_intercepts) {
    names_lambda <- c("Intercept")
    names_mu <- c("Intercept")
  }
  if (! rem_ran_ef) {
    if (length(attr(terms(nobars(formula_lambda)), "term.labels")) > 0) {
      names_lambda <- c(names_lambda, paste(attr(terms(nobars(formula_lambda)), "term.labels"), "lambda", sep = "_"))
    }
    if (length(attr(terms(nobars(formula_mu)), "term.labels")) > 0) {
      names_mu <- c(names_mu, paste(attr(terms(nobars(formula_mu)), "term.labels"), "mu", sep = "_"))
    }

    names(reduced_formulas_lambda) <- names_lambda
    names(reduced_formulas_mu) <- names_mu
  }

  if (type == 2) {
    print("Fitting full models")
    full_fits_lambda <- lapply(full_formulas_lambda, function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      return(fit_glmm(formula_tmp, data, mm))
    })
    full_fits_mu <- lapply(full_formulas_mu, function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      return(fit_glmm(formula_tmp, data, mm))
    })
    print("Fitting reduced models")
    full_formulas_lambda_char <- sapply(full_formulas_lambda, function(x) {as.character(x)[3]})
    reduced_fits_lambda <- lapply(reduced_formulas_lambda, function(formula_tmp) {
      if (as.character(formula_tmp)[3] %in% full_formulas_lambda_char) {
        print(paste("Copying ", formula_tmp, sep = ""))
        which_model_equal <- which(full_formulas_lambda_char == as.character(formula_tmp)[3])
        return(full_fits_lambda[[which_model_equal]])
      } else {
        print(paste("Fitting ", formula_tmp, sep = ""))
        return(fit_glmm(formula_tmp, data, mm))
      }
    })

    full_formulas_mu_char <- sapply(full_formulas_mu, function(x) {as.character(x)[3]})
    reduced_fits_mu <- lapply(reduced_formulas_mu, function(formula_tmp) {
      if (as.character(formula_tmp)[3] %in% full_formulas_mu_char) {
        print(paste("Copying ", formula_tmp, sep = ""))
        which_model_equal <- which(full_formulas_mu_char == as.character(formula_tmp)[3])
        print(which_model_equal)
        return(full_fits_mu[[which_model_equal]])
      } else {
        print(paste("Fitting ", formula_tmp, sep = ""))
        return(fit_glmm(formula_tmp, data, mm))
      }

    })
  } else {
    # fit reduced models
    print("Fitting reduced models")
    reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), function(formula_tmp) {
      print(paste("Fitting ", formula_tmp, sep = ""))
      return(fit_glmm(formula_tmp, data, mm))
    })
  }

  if (type == 2) return(list("reduced_fits_lambda" = reduced_fits_lambda,
                             "full_fits_lambda" = full_fits_lambda,
                             "reduced_fits_mu" = reduced_fits_mu,
                             "full_fits_mu" = full_fits_mu))
  else return(reduced_fits)
}


compute_LRTs <- function(fit_obj = NULL, formula_mu, formula_lambda, dv, data,
                         mm, type = 3, test_intercepts = F, test_ran_ef = F) {
  # only removes fixed effect, corresponding random slopes stay in the reduced model

  if (! type %in% c(2, 3)) stop("Please set type to 2 or 3. Returning NULL.")

  submodels <- fit_submodels(formula_mu, formula_lambda, dv, data, mm, type, test_intercepts,
                             rem_ran_ef = test_ran_ef)
  if(is.null(submodels)) return(NULL)
  if (type == 3) {
    reduced_fits <- submodels
    LL_full <- stats::logLik(fit_obj)
    LRT_results <- sapply(reduced_fits, function(fit_tmp) {

      LL_reduced <- stats::logLik(fit_tmp)
      chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
      df <- stats::df.residual(fit_tmp) - stats::df.residual(fit_obj)
      if (test_ran_ef) p_value <- pchisqmix(q = chisq, df = 1, lower.tail = F, mix = 0.5)
      else p_value <- pchisq(q = chisq, df = df, lower.tail = F)
      return(data.frame(
        # columns names inspired by afex
        "deviance_full" = -2 * LL_full,
        "deviance_reduced" = -2 * LL_reduced,
        "df.LRT" = df,
        "Chisq" = chisq,
        "p.value" = p_value
      ))
    })

    to_return <- list(
      "reduced_fits" = reduced_fits,
      "LRTs" = t(LRT_results)
    )
  }

  else if (type == 2) {
    #return(submodels)
    reduced_fits_lambda <- submodels[["reduced_fits_lambda"]]
    reduced_fits_mu <- submodels[["reduced_fits_mu"]]
    # full fit for the highest-order effects is the already-fitted full model
    # full_fits_lambda <- append(submodels[["full_fits_lambda"]], fit_obj)
    full_fits_lambda <- submodels[["full_fits_lambda"]]
    full_fits_lambda[[length(full_fits_lambda) + 1]] <- fit_obj
    full_fits_mu <- submodels[["full_fits_mu"]]
    full_fits_mu[[length(full_fits_mu) + 1]] <- fit_obj
    # assigns <- attr(mm[["lambda"]], "assign")

    # Compute Lambda LRTs if there are parameters to test
    if (length(reduced_fits_lambda) > 0) {
      orders_lambda <- attr(terms(lme4::nobars(formula_lambda)), "order")
      if (test_intercepts) orders_lambda <- c(0, orders_lambda)
      LRTs_lambda <- sapply(1:length(reduced_fits_lambda), function(fit_ind) {
        fit_tmp <- reduced_fits_lambda[[fit_ind]]
        LL_reduced <- stats::logLik(fit_tmp)
        LL_full <- stats::logLik(full_fits_lambda[[orders_lambda[fit_ind] + as.numeric(test_intercepts)]])

        chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
        p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
        df <- stats::df.residual(fit_tmp) - stats::df.residual(full_fits_lambda[[orders_lambda[fit_ind] + as.numeric(test_intercepts)]])
        return(data.frame(
          # columns names inspired by afex
          "deviance_full" = -2 * LL_full,
          "deviance_reduced" = -2 * LL_reduced,
          "df.LRT" = df,
          "Chisq" = chisq,
          "p.value" = p_value
        ))
      })
      colnames(LRTs_lambda) <- names(reduced_fits_lambda)
    } else LRTs_lambda <- NULL

    # Compute Mu LRTs if there are parameters to test
    if (length(reduced_fits_mu) > 0) {
      orders_mu <- attr(terms(lme4::nobars(formula_mu)), "order")
      if (test_intercepts) orders_mu <- c(0, orders_mu)
      LRTs_mu <- sapply(1:length(reduced_fits_mu), function(fit_ind) {
        fit_tmp <- reduced_fits_mu[[fit_ind]]
        LL_reduced <- stats::logLik(fit_tmp)
        LL_full <- stats::logLik(full_fits_mu[[orders_mu[fit_ind] + as.numeric(test_intercepts)]])
        chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
        p_value <- pchisq(q = chisq, df = 1, lower.tail = F)
        df <- stats::df.residual(fit_tmp) - stats::df.residual(full_fits_mu[[orders_mu[fit_ind] + as.numeric(test_intercepts)]])
        return(data.frame(
          # columns names inspired by afex
          "deviance_full" = -2 * LL_full,
          "deviance_reduced" = -2 * LL_reduced,
          "df.LRT" = df,
          "Chisq" = chisq,
          "p.value" = p_value
        ))
      })
      colnames(LRTs_mu) <- names(reduced_fits_mu)
    } else LRTs_mu <- NULL

    LRT_results <- cbind(LRTs_lambda, LRTs_mu)
    reduced_fits <- c(reduced_fits_lambda, reduced_fits_mu)
    full_fits <- c(full_fits_lambda, full_fits_mu)

    to_return <- list(
      "LRTs" = t(LRT_results),
      "reduced_fits" = reduced_fits,
      "full_fits" = full_fits
    )
  }
  # compute LRTs
  return(to_return)
}


pchisqmix <- function(q, df, mix, lower.tail = TRUE) {
  df_vec <- rep(df, length(q))
  mix_vec <- rep(mix, length(q))
  upper <- stats::pchisq(q = q, df = df, lower.tail = lower.tail)
  lower <- ifelse(df == 1, if (lower.tail) 1 else 0,
                  pchisq(q, df-1, lower.tail = lower.tail))
  return(mix * lower + (1 - mix) * upper)
}
