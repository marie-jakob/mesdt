fit_submodels <- function(formula_mu, formula_lambda, dv, data, mm, type = 3, test_intercepts = F) {
  # Default behavior: generate reduced models for all _variables_ (not factors) in the formula
  # -> correspond to multiple model parameters for factors with more than two levels

  if (test_intercepts) {
    #if (type == 2) {
    #  message("Cannot test intercepts for type II sums of squares (there's nothing left in the model).
    #                        Returning NULL.")
    #  return()
    #}
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
    #if (test_intercepts) orders <- c(0, orders)
    # no need to fit new full models if there are only main effects in the model
    if ((max(orders) < 2) & ! test_intercepts) full_formulas_lambda <- NULL
    else {
      if (test_intercepts) n_orders_lambda <- 0:(max(orders) - 1)
      else n_orders_lambda <- 1:(max(orders) - 1)
      full_formulas_lambda <- lapply(n_orders_lambda, function(x) {
        # add 0 for the intercept (which is dropped since it has no order)
        to_remove <- which(c(0, orders[assigns]) > x)
        if (length(to_remove) == dim(mm[["lambda"]])[2]) to_remove = Inf
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = F)
      })
    }
    # reduced_formulas are constructed for all predictors in the model for both mu and lambda
    # -> in order of the order of terms in the model matrix
    reduced_formulas_lambda <- lapply(range_lambda, function(x) {
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x))
      if (length(to_remove) == dim(mm[["lambda"]])[2]) to_remove = Inf
      print(paste("to_remove lambda: ", to_remove))
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              param_idc = to_remove, remove_from_mu = F)
    })

    assigns <- attr(mm[["mu"]], "assign")
    orders <- attr(terms(lme4::nobars(formula_mu)), "order")
    if ((max(orders) < 2) & ! test_intercepts) full_formulas_mu <- NULL
    else {
      if (test_intercepts) n_orders_mu <- 0:(max(orders) - 1)
      else n_orders_mu <- 1:(max(orders) - 1)

      full_formulas_mu <- lapply(n_orders_mu, function(x) {
        to_remove <- which(c(0, orders[assigns]) > x)
        print(paste("to_remove lambda: ", to_remove))
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                param_idc = to_remove, remove_from_mu = T)
      })
    }

    reduced_formulas_mu <- lapply(range_mu, function(x) {
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x))
      if (length(to_remove) == dim(mm[["mu"]])[2]) to_remove = Inf
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
    print(paste("length(full_fits_lambda): ", length(full_fits_lambda)))
    orders_lambda <- attr(terms(lme4::nobars(formula_lambda)), "order")
    if (test_intercepts) orders_lambda <- c(0, orders_lambda)
    LRTs_lambda <- sapply(1:length(reduced_fits_lambda), function(fit_ind) {
      print(fit_ind)
      print(orders_lambda[fit_ind] + as.numeric(test_intercepts))
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

    print(paste("length(full_fits_mu)", length(full_fits_mu)))
    print(paste("length(reduced_fits_mu)", length(reduced_fits_mu)))
    orders_mu <- attr(terms(lme4::nobars(formula_mu)), "order")
    if (test_intercepts) orders_mu <- c(0, orders_mu)
    LRTs_mu <- sapply(1:length(reduced_fits_mu), function(fit_ind) {
      fit_tmp <- reduced_fits_mu[[fit_ind]]
      LL_reduced <- stats::logLik(fit_tmp)
      LL_full <- stats::logLik(full_fits_mu[[orders_mu[fit_ind] + as.numeric(test_intercepts)]])
      print(paste("orders_lambda[fit_ind] + as.numeric(test_intercepts)", orders_lambda[fit_ind] + as.numeric(test_intercepts)))
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

    LRT_results <- cbind(LRTs_lambda, LRTs_mu)
  }


  # compute LRTs
  return(list(
    #"reduced_fits" = c(reduced_fits_lambda, reduced_fits_mu),
    "LRTs" = t(LRT_results)
  ))
}

