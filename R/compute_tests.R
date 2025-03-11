fit_submodels <- function(formula_mu, formula_lambda, dv, data, mm, type = 3, test_intercepts = F,
                          rem_ran_ef = F, correlate_sdt_params = T, cl = NULL,
                          test_params_mu = "all", test_params_lambda = "all",
                          control = NULL) {
  # Default behavior: generate reduced models for all _variables_ (not factors) in the formula
  # -> correspond to multiple model parameters for factors with more than two levels
  if (rem_ran_ef) {
    # range_lambda and range_mu are now list (to accommodate multiple random effects)
    range_lambda <- list(); range_mu <- list();
    for (i in 1:length(mm)) {
      if (grepl("rdm_lambda", names(mm)[i])) {
        if (test_intercepts) {
          range_lambda[[names(mm)[i]]] <- 0:max(attr(mm[[names(mm)[i]]], "assign"))
        } else {
          if (length(attr(mm[[names(mm)[i]]], "assign")) == 1) range_lambda[[names(mm)[i]]] <- NULL
          else range_lambda[[names(mm)[i]]] <- 1:max(attr(mm[[names(mm)[i]]], "assign"))
        }

      } else if (grepl("rdm_mu", names(mm)[i])) {
        if (test_intercepts) {
          range_mu[[names(mm)[i]]] <- 0:max(attr(mm[[names(mm)[i]]], "assign"))
        } else {
          if (length(attr(mm[[names(mm)[i]]], "assign")) == 1) range_mu[[names(mm)[i]]] <- NULL
          else range_mu[[names(mm)[i]]] <- 1:max(attr(mm[[names(mm)[i]]], "assign"))
        }
      }
    }
    if (length(range_lambda) == 0) range_lambda <- NULL
    if (length(range_mu) == 0) range_mu <- NULL
  } else {


    if (test_intercepts) {
      range_lambda <- list(0:max(attr(mm[["lambda"]], "assign")))
      range_mu <- list(0:max(attr(mm[["mu"]], "assign")))
      print(range_lambda)
      print(range_mu)
    } else {
      if (length(attr(mm[["lambda"]], "assign")) == 1) range_lambda <- c()
      else range_lambda <- list(1:max(attr(mm[["lambda"]], "assign")))
      if (length(attr(mm[["mu"]], "assign")) == 1) range_mu <- c()
      else range_mu <- list(1:max(attr(mm[["mu"]], "assign")))
    }


    # Update range of to-be-tested parameters if only some parameters are to be tested
    if (is.null(test_params_mu)) {
      labels_mu <- NULL; test_mu_labels <- NULL; range_mu <- NULL;
    } else if (test_params_mu != "all") {
      labels_mu <- attr(stats::terms.formula(lme4::nobars(formula_mu)), "term.labels")
      test_mu_labels <- attr(stats::terms.formula(test_params_mu), "term.labels")
      if (! all_terms_in(test_mu_labels, labels_mu)) stop("Only parameters that are in the model can be tested.")
      range_mu <- ifelse(test_intercepts,
                         list(range_mu[[1]][c(1, which_terms_in(reference_terms = labels_mu,
                                                                terms_to_check = test_mu_labels) + 1)]),
                         list(range_mu[[1]][which_terms_in(reference_terms = labels_mu,
                                                           terms_to_check = test_mu_labels)]))
    }
    if (is.null(test_params_lambda)) {
      labels_lambda <- NULL; test_lambda_labels <- NULL; range_lambda <- NULL;
    } else if (test_params_lambda != "all") {
      labels_lambda <- attr(stats::terms.formula(lme4::nobars(formula_lambda)), "term.labels")
      test_lambda_labels <- attr(stats::terms.formula(test_params_lambda), "term.labels")
      if (! all_terms_in(test_lambda_labels, labels_lambda)) stop("Only parameters that are in the model can be tested.")
      range_lambda <- ifelse(test_intercepts,
                             list(range_lambda[[1]][c(1, which_terms_in(reference_terms = labels_lambda,
                                                                        terms_to_check = test_lambda_labels) + 1)]),
                             list(range_lambda[[1]][which_terms_in(reference_terms = labels_lambda,
                                                                   terms_to_check = test_lambda_labels)]))
    }
    print(range_lambda)
    print(range_mu)

    if (length(range_lambda) > 0) {
      names(range_lambda) <- "lambda"
    }
    if (length(range_mu) > 0) {
      names(range_mu) <- "mu"
    }
    if (length(range_lambda) == 0 & length(range_mu) == 0) {
      message("Nothing to test in the model. Returning NULL.")
      return(NULL)
    }
  }

  if (type == 3) {
    reduced_formulas_lambda <- list()
    for (i in 1:length(range_lambda)) {
      forms_tmp <- lapply(range_lambda[[i]], function(x) {
        name_tmp <- names(range_lambda)[i]
        to_remove <- list()
        to_remove[[name_tmp]] <- as.numeric(which(attr(mm[[name_tmp]], "assign") == x))
        if (length(to_remove[[name_tmp]]) == dim(mm[[name_tmp]])[2]) to_remove[[name_tmp]] = Inf
        return(construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                      to_remove = to_remove, correlate_sdt_params = correlate_sdt_params))
      })
      reduced_formulas_lambda <- append(reduced_formulas_lambda, forms_tmp)
    }
    reduced_formulas_mu <- list()
    for (i in 1:length(range_mu)) {
      forms_tmp <- lapply(range_mu[[i]], function(x) {
        name_tmp <- names(range_mu)[i]
        to_remove <- list()
        to_remove[[name_tmp]] <- as.numeric(which(attr(mm[[name_tmp]], "assign") == x))
        if (length(to_remove[[name_tmp]]) == dim(mm[[name_tmp]])[2]) to_remove[[name_tmp]] = Inf
        return(construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                       to_remove = to_remove, correlate_sdt_params = correlate_sdt_params))
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

    if (is.null(test_params_lambda)) full_formulas_lambda <- NULL
    if ((max(c(0, orders)) + as.numeric(test_intercepts)) < 2) full_formulas_lambda <- NULL
    else {
      if (test_intercepts) n_orders_lambda <- 0:(max(orders) - 1)
      else n_orders_lambda <- 1:(max(orders) - 1)
      full_formulas_lambda <- lapply(n_orders_lambda, function(x) {
        # add 0 for the intercept (which is dropped since it has no order)
        to_remove <- list()
        to_remove[["lambda"]] <- as.numeric(which(c(0, orders[assigns]) > x))
        if (length(to_remove[["lambda"]]) == dim(mm[["lambda"]])[2]) to_remove[["lambda"]] <- Inf
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                to_remove = to_remove, correlate_sdt_params = correlate_sdt_params)
      })
    }
    # reduced_formulas are constructed for all predictors in the model for both mu and lambda
    # -> in order of the order of terms in the model matrix
    reduced_formulas_lambda <- lapply(range_lambda[[1]], function(x) {
      #print(x)
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- list()
      to_remove[["lambda"]] <- as.numeric(c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x)))
      if (length(to_remove[["lambda"]]) == dim(mm[["lambda"]])[2]) to_remove[["lambda"]] <- Inf
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              to_remove = to_remove, correlate_sdt_params = correlate_sdt_params)
    })

    assigns <- attr(mm[["mu"]], "assign")
    orders <- attr(terms(lme4::nobars(formula_mu)), "order")

    if (is.null(test_params_mu)) full_formulas_mu <- NULL
    if ((max(c(0, orders)) + as.numeric(test_intercepts)) < 2) full_formulas_mu <- NULL
    else {
      if (test_intercepts) n_orders_mu <- 0:(max(orders) - 1)
      else n_orders_mu <- 1:(max(orders) - 1)

      full_formulas_mu <- lapply(n_orders_mu, function(x) {
        to_remove <- list()
        to_remove[["mu"]] <- as.numeric(which(c(0, orders[assigns]) > x))
        if (length(to_remove[["mu"]]) == dim(mm[["mu"]])[2]) to_remove[["mu"]] <- Inf
        construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                                to_remove = to_remove, correlate_sdt_params = correlate_sdt_params)
      })
    }
    reduced_formulas_mu <- lapply(range_mu[[1]], function(x) {
      if (x == 0) order_of_x <- 0
      else order_of_x <- orders[assigns][x]
      to_remove <- list()
      to_remove[["mu"]] <- as.numeric(c(which(c(0, orders[assigns]) > order_of_x), which(assigns == x)))
      if (length(to_remove[["mu"]]) == dim(mm[["mu"]])[2]) to_remove[["mu"]] <- Inf
      construct_glmer_formula(formula_mu, formula_lambda, dv = dv, mm = mm,
                              to_remove = to_remove, correlate_sdt_params = correlate_sdt_params)
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
    if (length(attr(terms(lme4::nobars(formula_lambda)), "term.labels")) > 0) {
      if (! is.null(test_params_lambda)) {
        if (test_params_lambda != "all") names_lambda <- paste(c(names_lambda, test_lambda_labels), "lambda", sep = "_")
        else names_lambda <- c(names_lambda, paste(attr(terms(lme4::nobars(formula_lambda)), "term.labels"), "lambda", sep = "_"))
      }
    }

    if (! is.null(test_params_mu)) {
      if (length(attr(terms(lme4::nobars(formula_mu)), "term.labels")) > 0) {
        if (test_params_mu != "all") names_mu <- paste(c(names_mu, test_mu_labels), "mu", sep = "_")
        else names_mu <- c(names_mu, paste(attr(terms(lme4::nobars(formula_mu)), "term.labels"), "mu", sep = "_"))
      }
    }
    if (length(reduced_formulas_lambda) != 0) names(reduced_formulas_lambda) <- names_lambda
    if (length(reduced_formulas_mu) != 0) names(reduced_formulas_mu) <- names_mu
  }

  if (type == 2) {
    if (is.null(cl)) {
      print("No cluster")

      full_fits <- lapply(c(full_formulas_lambda, full_formulas_mu),
                          fit_glmm, data = data, mm = mm, control = control)
      # for reduced fits: check if the fit is already in the full_fits

      full_formulas_char <- sapply(c(full_formulas_lambda, full_formulas_mu), function(x) {as.character(x)[3]})

      reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), function(formula_tmp) {
        if (as.character(formula_tmp)[3] %in% full_formulas_char) {
          #print(paste("Copying ", formula_tmp, sep = ""))
          which_model_equal <- which(full_formulas_char == as.character(formula_tmp)[3])
          #print(which_model_equal)
          return(full_fits[[which_model_equal]])
        } else return(fit_glmm(formula_tmp, data, mm, control))
      })
      all_fits <- c(full_fits, reduced_fits)
    } else {
      print("cluster")
      throwaway <- parallel::clusterExport(cl = cl,
                                           varlist = c("data", "mm", "fit_glmm"),
                                           env = environment())
      all_fits <- unlist(parallel::clusterApplyLB(cl,
                                                  c(full_formulas_lambda,
                                                    full_formulas_mu,
                                                    reduced_formulas_lambda,
                                                    reduced_formulas_mu),
                                                    fit_glmm,
                                                    data = data, mm = mm, control = control))
      names(all_fits) <- names(c(full_formulas_lambda,
                                 full_formulas_mu,
                                 reduced_formulas_lambda,
                                 reduced_formulas_mu))

    }
    if (!length(full_formulas_lambda) == 0) {
      full_fits_lambda <- all_fits[1:length(full_formulas_lambda)]
    } else full_fits_lambda <- NULL
    if (!length(full_formulas_mu) == 0) {
      full_fits_mu <- all_fits[(length(full_formulas_lambda) + 1) : (length(full_formulas_lambda) + length(full_formulas_mu))]
    } else full_fits_mu <- NULL
    if (!length(reduced_formulas_lambda) == 0) {
      reduced_fits_lambda <- all_fits[(length(full_formulas_lambda) + length(full_formulas_mu) + 1):
                                        (length(all_fits) - length(reduced_formulas_mu))]
    } else reduced_fits_lambda <- NULL
    if (!length(reduced_formulas_mu) == 0) {
      reduced_fits_mu <- all_fits[(length(all_fits) - length(reduced_formulas_mu) + 1) : length(all_fits)]
    } else reduced_fits_mu <- NULL

    return(list(
      "full_fits_lambda" = full_fits_lambda,
      "full_fits_mu" = full_fits_mu,
      "reduced_fits_lambda" = reduced_fits_lambda,
      "reduced_fits_mu" = reduced_fits_mu
    ))
  }

  else {
    if (is.null(cl)) {
      print("No cluster")
      reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), fit_glmm,
                             data = data, mm = mm, control = control)
      print("reduced_fits done")
    } else {
      print("cluster")
      throwaway <- parallel::clusterExport(cl = cl,
                                           varlist = c("data", "mm", "fit_glmm", "control"),
                                           env = environment())
      reduced_fits <- unlist(parallel::clusterApplyLB(cl,
                                                      c(reduced_formulas_lambda, reduced_formulas_mu),
                                                      fit_glmm,
                                                      data = data, mm = mm, control = control))
      names(reduced_fits) <- names(c(reduced_formulas_lambda, reduced_formulas_mu))
    }
    return(reduced_fits)
  }
}

#' Compute Likelihood Ratio Tests or tests based on parametric bootstrapping
#' for a fitted SDT model
#'
#' @param fit_obj An lme4 or glmmTMB fit object containing the full fit and
#'  coefficients that should be tested
#' @param formula_mu the corresponding formula for sensitivity
#' @param formula_lambda the corresponding formula for response bias
#' @param dv the name of the dependent variable in the data
#' @param data a data frame
#' @param trial_type_var variable in the dataframe coding the type of trial (signal vs. noise)
#' @param tests type of tests that should be computed ("LRT" -> likelihood ratio tests,
#' "bootstrap" = parametric bootstrap)
#' @param nsim number of simulated datasets for bootstrapping
#' @param mm model matrices (optional)
#' @param type type of tests (only relevant for likelihood ratio tests and
#'  parametric bootstrapping)
#' @param test_intercepts boolean indicating if intercepts for sensitivity and
#'  response bias should be tested
#' @param test_ran_ef boolean indicating whether random (T) or fixed effects
#'  should be tested (F)
#' @param correlate_sdt_params boolean indicating if correlations between SDT
#'  parameters should be estimated
#' @param test_params_mu which coefficients on sensitivity should be tested?
#'  (default is "all")
#' @param test_params_lambda which coefficients on response bias should be tested?
#'  (default is "all")
#'
#' @export
compute_tests <- function(mlsdt_fit, data,
                          cl = NULL, tests = "LRT", nsim = 1000,
                          mm = NULL, type = 3, test_intercepts = F, test_ran_ef = F,
                          correlate_sdt_params = T, test_params_mu = "all",
                          test_params_lambda = "all",
                          control = NULL) {

  fit_obj <- mlsdt_fit$fit_obj
  formula_mu <- mlsdt_fit$formula_mu
  formula_lambda <- mlsdt_fit$formula_lambda
  dv <- mlsdt_fit$dv
  trial_type_var <- mlsdt_fit$trial_type_var

  # only removes fixed effect, corresponding random slopes stay in the reduced model

  if (is.null(mm)) {
    mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var = trial_type_var)
  }
  if (test_ran_ef & type != 3) {
    stop("Only type III sums of squares are available for testing random effects.")
  }

  if (! type %in% c(2, 3)) stop("Please set type to 2 or 3. Returning NULL.")
  if (! tests %in% c("LRT", "bootstrap")) stop('Only likelihood ratio tests (type = "LRT") and parametric bootstrapping(type = "bootstrap")
                                         are implemented')

  #if (is.null(seed)) seed <- sample(1:1e6, 1)
  # TODO: test compatibility of input arguments

  submodels <- fit_submodels(formula_mu, formula_lambda, dv, data, mm, type, test_intercepts = test_intercepts,
                             rem_ran_ef = test_ran_ef, correlate_sdt_params = correlate_sdt_params,
                             cl = cl, control = control,
                             test_params_mu = test_params_mu, test_params_lambda = test_params_lambda)
  if(is.null(submodels)) return(NULL)
  if (type == 3) {
    reduced_fits <- submodels
    LRT_results <- sapply(reduced_fits, compute_lrt, fit_full = fit_obj, test_ran_ef = test_ran_ef)

    if (tests == "bootstrap") {
      pb_objects <- lapply(reduced_fits, function(fit_tmp) {
        boot_tmp <- compute_parametric_bootstrap_test(fit_obj, fit_tmp, dv = dv, data = data,
                                                      mm = mm, nsim = nsim, cl = cl, control = control)
        return(boot_tmp)
      })
      boot_table <- sapply(pb_objects, function(x) { return(data.frame(x$test)[2, ]) })
      to_return <- list(
        "reduced_fits" = reduced_fits,
        "LRTs" = t(LRT_results),
        "pb_tests" = t(boot_table),
        "pb_objects" = pb_objects
      )
    } else {
      to_return <- list(
        "reduced_fits" = reduced_fits,
        "LRTs" = t(LRT_results)
      )
    }
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
      #print(reduced_fits_lambda)
      #print(lme4::nobars(formula_lambda))
      orders_lambda <- attr(terms(lme4::nobars(formula_lambda)), "order")
      if (test_intercepts) orders_lambda <- c(0, orders_lambda)
      LRTs_lambda <- sapply(1:length(reduced_fits_lambda), function(fit_ind) {
        fit_tmp <- reduced_fits_lambda[[fit_ind]]
        LL_reduced <- stats::logLik(fit_tmp)
        LL_full <- stats::logLik(full_fits_lambda[[orders_lambda[fit_ind] + as.numeric(test_intercepts)]])
        chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
        p_value <- stats::pchisq(q = chisq, df = 1, lower.tail = F)
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
      if (tests == "bootstrap") {
        pb_objects_lambda <- lapply(1:length(reduced_fits_lambda), function(fit_ind) {
          fit_tmp <- reduced_fits_lambda[[fit_ind]]
          fit_full <- full_fits_lambda[[orders_lambda[fit_ind] + as.numeric(test_intercepts)]]
          boot_tmp <- compute_parametric_bootstrap_test(fit_full, fit_tmp, data = data,
                                                        mm = mm, dv = dv, nsim = nsim, cl = cl, control = control)
          return(boot_tmp)
        })
        boot_table_lambda <- sapply(pb_objects_lambda, function(x) { return(data.frame(x$test)[2, ]) })
        colnames(boot_table_lambda) <- names(reduced_fits_lambda)
      }
    } else {
      LRTs_lambda <- NULL
      if (tests == "bootstrap") boot_table_lambda <- NULL; pb_objects_lambda <- NULL;
    }

    # Compute Mu LRTs if there are parameters to test
    if (length(reduced_fits_mu) > 0) {
      orders_mu <- attr(terms(lme4::nobars(formula_mu)), "order")
      if (test_intercepts) orders_mu <- c(0, orders_mu)
      LRTs_mu <- sapply(1:length(reduced_fits_mu), function(fit_ind) {
        fit_tmp <- reduced_fits_mu[[fit_ind]]
        LL_reduced <- stats::logLik(fit_tmp)
        LL_full <- stats::logLik(full_fits_mu[[orders_mu[fit_ind] + as.numeric(test_intercepts)]])
        chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
        p_value <- stats::pchisq(q = chisq, df = 1, lower.tail = F)
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
      if (tests == "bootstrap") {
        pb_objects_mu <- lapply(1:length(reduced_fits_mu), function(fit_ind) {
          fit_tmp <- reduced_fits_mu[[fit_ind]]
          fit_full <- full_fits_mu[[orders_mu[fit_ind] + as.numeric(test_intercepts)]]
          boot_tmp <- compute_parametric_bootstrap_test(fit_full, fit_tmp, data = data,
                                                        mm = mm, dv = dv, nsim = nsim, cl = cl, control = control)
          return(boot_tmp)
        })
        boot_table_mu <- sapply(pb_objects_mu, function(x) { return(data.frame(x$test)[2, ]) })
        colnames(boot_table_mu) <- names(reduced_fits_mu)
      }

    } else {
      LRTs_mu <- NULL
      if (tests == "bootstrap") boot_table_mu <- NULL; pb_objects_mu <- NULL;
    }

    LRT_results <- cbind(LRTs_lambda, LRTs_mu)
    reduced_fits <- c(reduced_fits_lambda, reduced_fits_mu)
    full_fits <- c(full_fits_lambda, full_fits_mu)

    if (tests == "bootstrap") {

      pb_objects <- c(pb_objects_lambda, pb_objects_mu)
      boot_table <- cbind(boot_table_lambda, boot_table_mu)
      to_return <- list(
        "LRTs" = t(LRT_results),
        "reduced_fits" = reduced_fits,
        "full_fits" = full_fits,
        "pb_tests" = t(boot_table),
        "pb_objects" = pb_objects
      )
    } else {
      to_return <- list(
        "LRTs" = t(LRT_results),
        "reduced_fits" = reduced_fits,
        "full_fits" = full_fits
      )
    }
  }
  return(to_return)
}


# Compute LRT of one fixed or random effect based on the given full model
# and submodel fit
compute_lrt <- function(fit_red, fit_full, test_ran_ef) {
  LL_full <- stats::logLik(fit_full)
  LL_reduced <- stats::logLik(fit_red)
  chisq <- (-2) * (as.numeric(LL_reduced) - as.numeric(LL_full))
  df <- stats::df.residual(fit_red) - stats::df.residual(fit_full)
  if (test_ran_ef) p_value <- pchisqmix(q = chisq, df = df, lower.tail = F, mix = 0.5)
  else p_value <- stats::pchisq(q = chisq, df = df, lower.tail = F)
  return(data.frame(
    # columns names inspired by afex
    "deviance_full" = -2 * LL_full,
    "deviance_reduced" = -2 * LL_reduced,
    "df.LRT" = df,
    "Chisq" = chisq,
    "p.value" = p_value
  ))
}


# Built after pbkrtest code
# https://hojsgaard.github.io/pbkrtest/index.html
# supports lme4 and glmmTMB objects
# TODO: export?
compute_parametric_bootstrap_test <- function(large_model, small_model, data, mm, dv,
                                              nsim = 1000, cl = NULL, control = NULL) {
  #if (is.null(seed)) seed <- sample(1:1e6, 1)

  do_pb <- function(x) {
    dat_tmp <- data
    sim_dat_tmp <- stats::simulate(small_model)[[1]]
    # Do refitting manually
    dat_tmp[[dv]] <- sim_dat_tmp
    sim_fit_full <- fit_glmm(stats::formula(large_model), dat_tmp, mm, control)
    sim_fit_red <- fit_glmm(stats::formula(small_model), dat_tmp, mm, control)
    return(-2 * (stats::logLik(sim_fit_red) - stats::logLik(sim_fit_full)))
  }

  if (is.null(cl)) {
    parallel <- F
    LRs_boot <- sapply(1:nsim, do_pb)
  } else {
    parallel <- T
    # load variables on cluster
    throwaway <- parallel::clusterExport(cl = cl,
                            varlist = c("data", "mm", "fit_glmm", "dv", "small_model", "large_model", "control"),
                            env = environment())
    if (options("mlsdt.backend") == "glmmTMB") throwaway <- parallel::clusterCall(cl = cl, "require", package = "glmmTMB", character.only = T)
    # set seed on cluster
    LRs_boot <- unlist(parallel::clusterApplyLB(cl, 1:nsim, do_pb))
  }
  LR_emp <- -2 * (stats::logLik(small_model) - stats::logLik(large_model))
  df_lrt <- stats::df.residual(small_model) - stats::df.residual(large_model)

  # only use samples where the LR is positive (as in pbkrtest)
  LRs_boot <- LRs_boot[which(LRs_boot > 0)]
  p_boot <- (length(which(LRs_boot > LR_emp)) + 1) / (length(LRs_boot) + 1)

  test <- data.frame("stat" = rep(LR_emp, 2),
                     "df" = c(df_lrt, NA),
                     "p.value" = c(stats::pchisq(LR_emp, df = df_lrt, lower.tail = F), p_boot))
  rownames(test) <- c("LRT", "PBtest")
  return(list("test" = test,
              "requested_samples" = nsim,
              "used_samples" = length(LRs_boot),
              "ref" = LRs_boot,
              "n.extreme" = length(which(LRs_boot > LR_emp)),
              "parallel" = parallel))
}


pchisqmix <- function(q, df, mix, lower.tail = TRUE) {
  df_vec <- rep(df, length(q))
  mix_vec <- rep(mix, length(q))
  upper <- stats::pchisq(q = q, df = df, lower.tail = lower.tail)
  lower <- ifelse(df_vec == 1, if (lower.tail) 1 else 0,
                  stats::pchisq(q, df-1, lower.tail = lower.tail))
  return(mix_vec * lower + (1 - mix_vec) * upper)
}
