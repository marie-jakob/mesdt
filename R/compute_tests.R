#' Compute likelihood ratio tests and parametric bootstrapping tests for a
#' (mixed-effects) SDT model.
#'
#' @param fit_obj `mesdt_fit` object containing the full fit and the
#'  coefficients that should be tested.
#' @param tests `character` specifying the type of tests that should be
#' computed ("LRT" for likelihood ratio tests, "bootstrap" for parametric
#' bootstrapping tests, see Details).
#' @param type `integer` or `character` specifying the type of tests
#'  (default 3 / "III", 2 / "II") that are computed (see Details).
#' @param test_intercepts `boolean` indicating if intercepts for
#' discriminability and response bias should be tested.
#' @param test_discriminability `formula` or `character` indicating which
#'  coefficients on discriminability should be tested (default is "all").
#'  See Details.
#' @param test_response_bias `formula` or `character` indicating which
#'  coefficients on response bias should be tested (default is "all").
#'  See Details.
#' @param nsim `integer` specifying the number of simulated datasets for
#'  bootstrapping.
#' @param cl name of a cluster to distribute computation across multiple cores
#' @param control list containing optional control arguments that are included
#'  in the `glmmTMB()` or `glmer()` call (see Details).

#'
#' @return An object of class `mesdt_test`, containing the fitted model
#' (`$fit_obj`), fit objects for the reduced fits (`$reduced_fits`), fit objects
#' for the full fits (`$full_fits`, only for type II tests), a dataframe
#' containing results of likelihood ratio tests (`$LRT_results`), the type of
#' tests that were computed (`$type`). For parametric bootstrap tests, there are
#' additional slots for the results of those tests (`$pb_test_results`), further
#' information about the tests (`$pb_objects`, similar to the structure
#' returned by the `pbkrtest` package), and the random seed that was used to
#' simulate the data (`$seed`).
#'
#' @description
#' Computes likelihood ratio tests and parametric bootstrapping tests for a
#' (mixed-effects) SDT model. The function allows the user to either test all
#' effects on sensitivity, response bias or both
#' (`tests_discriminability = "all"`, `tests_response_bias = "all"`,
#' the default), to specify
#' effects to be tested (by providing formulas, see Examples), or to don't
#' test any effect on one of the SDT parameters
#' (`tests_discriminability = NULL` or `tests_response_bias = NULL`).
#' Intercepts for response bias and discriminability can optionally be tested as
#' well by setting `test_intercepts = TRUE`. A `print` method is provided,
#' showing the test results; Everything else can be directly extracted
#' from the returned `mesdt_test` object.
#'
#' Both methods compute tests for the given effects by comparing a full model,
#' containing the fixed effect of interest with a restricted model, where the
#' fixed effect of interest is fixed to zero. For likelihood ratio tests
# ("lrt", the default), the resulting likelihood ratio is taken as an empirical
#' Chi^2 value and p-values are computed with a Chi^2 test. For parametric
#' bootstrapping, p-values are computed based on a simulated reference
#' distribution of the restricted model (see Details). The function implements
#' type III (the default) and type II tests, as specified with the `type`
#' argument.
#'
#' @details
#'
#' ## Test Logic
#'  Likelihood ratio tests (LRTs) and parametric bootstrapping tests (PBTs)
#'  are based on comparisons of nested models: A full model, containing the
#'  effect of interest, and a restricted model, where the fixed effect of
#'  interest is set to zero. The resulting difference in the deviance of both
#'  models is a likelihood ratio, which asymptotically follows a Chi^2
#'  distribution, with degrees of freedom equal to the difference in the
#'  numbers of parameters of both models. LRTs thus use this Chi^2 distribution
#'  as the reference distribution. However, because this is an _asymmetric_
#'  property of the likelihood ratio, it is usually only recommended for models
#'  and data with a sufficiently high number of levels for the random-effects
#'  grouping factors (e.g., > 50, as recommended in the afex package). For
#'  other cases, PBTs can be used (note, however, that the computation is
#'  computationally intensive): PBTs generate the reference distribution by
#'  repeatedly simulating data from the restricted model and estimating the
#'  full model on the simulated data.
#'
#'
#' ## Model Specification of Type II and III tests
#'  Type III and type II tests follow a similar logic as the types of sums of
#'  squares in analyses of variance: For type III (the default) the full model
#'  is always the given model and the restricted models are specified by
#'  only removing the fixed effect of interest from the model. In contrast,
#'  type II sums of squares take into account the order of the to-be-tested
#'  effect and aim to preserve the principle of marginality (i.e., models
#'  including higher-order effects such as interactions should always include
#'  all effects that are marginal to them). Thus, type II tests are obtained by
#'  comparing a model in which the to-be-tested effect and all higher-order
#'  effects are removed with a model in which only effects up to the effect of
#'  the to-be-tested parameter are included (similar to the logic applied in
#'  afex). Importantly, this logic is applied only to the model structure for
#'  the SDT parameter for which the effect is to be tested (e.g., when testing
#'  a main effect of a variable on response bias, all higher-order effects are
#'  removed from response bias, while the model for discriminability remains
#'  unchanged). Thus, marginality is preserved _within_ the respective SDT
#'  parameter, but higher-order interactions of to-be-tested effects can
#'  (if specified) be included on the respective other SDT parameter.
#'
#' ## Random effects of Restricted Models
#'  The hypothesis tests implemented here only concern fixed effects, meaning
#'  that the random effects structure remains unchanged. Thus, when testing
#'  effects of a variable x, the random slope for x remains in the restricted
#'  model, if it is present in the full model. This corresponds to what has
#'  been referred to as a _balanced null_ model (unlike a _strict null_ model,
#'  in which a possible corresponding random slope is removed from the
#'  restricted model as well). Substantively, these types of model comparison
#'  test whether the _mean population effect_ is significantly different from
#'  zero, regardless of variability in this effect according to random grouping
#'  factors. Thus, in case of a non-significant effect, individuals might still
#'  exhibit non-zero effects.
#'
#'
#' ## Computational Details
#'
#'  The `control` argument allows to pass all control arguments taken by
#'  `lme4` or `glmmTMB` (depending on the backend used; see the relevant
#'  documentation for details).
#'
#' @examples
#' \dontrun{
#' mod <- fit_mesdt(
#' discriminability ~ committee * emp_gender + (1 | id),
#' bias ~ committee * emp_gender + (committee | id),
#' data = debi3_sub,
#' trial_type = "status",
#' dv = "assessment"
#' )
#' # LRTs for all parameters
#' tests_all <- compute_tests(mod)
#' tests_all
#' # Tests intercepts as well
#' tests_all <- compute_tests(mod, test_intercepts = TRUE)
#' tests_all
#' # Only test effect of committee on response bias
#' tests_committee <- compute_tests(mod,
#'                                  tests_discriminability = NULL,
#'                                  tests_response_bias = ~ committee)
#' tests_committee
#' # parametric bootstrapping tests with a given cluster
#' require(parallel)
#' cl <- makeCluster(4, type = "SOCK")
#' tests_boot <- compute_tests(fit,
#'                             tests_response_bias = ~ committee,
#'                             tests_discriminability = NULL,
#'                             tests = "bootstrap",
#'                             nsim = 1000,
#'                             cl = cl,
#'                             seed = 51)
#' stopCluster(cl)
#' tests_boot
#'}
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterCall
#' @importFrom parallel clusterApplyLB
#' @importFrom parallel clusterExport
#' @export
compute_tests <- function(mesdt_fit, tests = "lrt",
                          type = 3, test_intercepts = F,
                          tests_discriminability = "all",
                          tests_response_bias = "all",
                          nsim = 1000, cl = NULL, control = NULL, seed = NULL) {
  ##### Check Input
  # TODO: give data as argument or take the data from the fitted model?
  # -> could lead to issues with extreme value
  test_ran_ef <- F

  if (! inherits(mesdt_fit, "mesdt_fit")) stop("Input 'mesdt_fit' must be of class 'mesdt_fit'.")
  type <- standardize_type_input(type)
  tests <- standardize_tests_input(tests)
  if (is.null(tests)) stop("Input 'tests' should be 'lrt' (for likelihood ratio tests) or 'bootstrap' (for tests based on parametric bootstrapping")
  if (! is.numeric(nsim)) stop("Input 'nsim' must be numeric.")
  if (round(nsim, 0) != nsim) stop("Input 'nsim' must be an integer.")
  if (is.null(type)) stop("Input 'type' should be 2 or 3.")

  if (! is.logical(test_intercepts)) stop("Input 'test_intercepts' must be of type 'logical'.")
  if (! is.logical(test_ran_ef)) stop("Input 'test_ran_ef' must be of type 'logical'.")

  if (! is.null(tests_discriminability)) {
    if (! typeof(tests_discriminability) == "language" & ! tests_discriminability == "all") stop("Input 'tests_discriminability' must be a formula or 'all'.")
    if (typeof(tests_discriminability) == "language") {
      # check that the formula does not contain random effects
      if (lme4::nobars(tests_discriminability) != tests_discriminability) stop("Input 'tests_discriminability' must not contain random effects.")
      else {
        terms_mod <- attr(terms(lme4::nobars(mesdt_fit$user_input$discriminability)), "term.labels")
        terms_input <- attr(terms(lme4::nobars(tests_discriminability)), "term.labels")
        if (! all(terms_input %in% terms_mod)) stop("Input 'tests_discriminability' contains terms that are not present in the fitted model. Please check your formula again.")
      }
    }
  }

  if (! is.null(tests_response_bias)) {
    if (! typeof(tests_response_bias) == "language" & ! tests_response_bias == "all") stop("Input 'tests_response_bias' must be a formula or 'all'.")
    if (typeof(tests_response_bias) == "language") {
      # check that the formula does not contain random effects
      if (lme4::nobars(tests_response_bias) != tests_response_bias) stop("Input 'tests_response_bias' must not contain random effects.")
      else {
        terms_mod <- attr(terms(lme4::nobars(mesdt_fit$user_input$bias)), "term.labels")
        terms_input <- attr(terms(lme4::nobars(tests_response_bias)), "term.labels")
        if (! all(terms_input %in% terms_mod)) stop("Input 'tests_response_bias' contains terms that are not present in the fitted model. Please check your formula again.")
      }
    }
  }

  fit_obj <- mesdt_fit$fit_obj
  data <- mesdt_fit$internal$data
  formula_mu <- mesdt_fit$user_input$discriminability
  formula_lambda <- mesdt_fit$user_input$bias
  dv <- mesdt_fit$user_input$dv
  trial_type_var <- mesdt_fit$user_input$trial_type_var
  correlate_sdt_params <- mesdt_fit$user_input$correlate_sdt_params
  distribution <- mesdt_fit$user_input$distribution

  mm <- mesdt_fit$internal$mm
  # only removes fixed effect, corresponding random slopes stay in the reduced model


  #if (is.null(mm)) {
  #  mm <- construct_modelmatrices(formula_mu, formula_lambda, data, trial_type_var = trial_type_var, distribution = distribution)[["mm"]]
  #}
  if (test_ran_ef & type != 3) {
    stop("Only type III sums of squares are available for testing random effects.")
  }

  if (! type %in% c(2, 3)) stop("Please set type to 2 or 3. Returning NULL.")

  # check if backend is the same as for the fitted model
  #print(paste("backend:", mesdt_fit$user_input$backend))
  if (is.null(cl)) {
    # if not, set the correct backend and notify the user
    if (options("mesdt.backend") != mesdt_fit$user_input$backend) {
      message(paste("Model was fitted using", mesdt_fit$user_input$backend, "but the current
                  backend is ", options("mesdt.backend"), ". Setting mesdt.backend
                  to ", mesdt_fit$user_input$backend))
      options("mesdt.backend" = mesdt_fit$user_input$backend)
    }
    # same when a cluster is used
  } else {
    backend_cl <- unname(unlist(parallel::clusterEvalQ(cl, options("mesdt.backend"))))[1]
    if (is.null(backend_cl)) {
      message(paste("No backend was set for the cluster. Setting mesdt.backend on the cluster and locally to", mesdt_fit$user_input$backend, "which was used to fit the supplied model."))
      throwaway <- parallel::clusterCall(cl, function() { options("mesdt.backend" = mesdt_fit$user_input$backend) } )
      # print(paste("mesdt.backend: ", mesdt_fit$user_input$backend))
      options("mesdt.backend" = mesdt_fit$user_input$backend)
    } else {
      if (backend_cl != mesdt_fit$user_input$backend) {
        message(paste("Model was fitted using", mesdt_fit$user_input$backend, "but the current
                  backend on the supplied cluster is ", backend_cl,
                  ". Setting msldt.backend on the cluster and locally to ", mesdt_fit$user_input$backend))
        throwaway <- parallel::clusterCall(cl, function() { options("mesdt.backend" = mesdt_fit$user_input$backend) } )
        options("mesdt.backend" = mesdt_fit$user_input$backend)
      }
    }

  }

  if (is.null(seed)) seed <- sample(1:1e8, 1)
  # TODO: test compatibility of input arguments

  submodels <- fit_submodels(formula_mu, formula_lambda, dv, data, mm, type, test_intercepts = test_intercepts,
                             rem_ran_ef = test_ran_ef, correlate_sdt_params = correlate_sdt_params,
                             distribution = distribution,
                             cl = cl, control = control,
                             tests_discriminability = tests_discriminability, tests_response_bias = tests_response_bias)
  if(is.null(submodels)) return(NULL)
  if (type == 3) {
    reduced_fits <- submodels
    LRT_results <- sapply(reduced_fits, compute_lrt, fit_full = fit_obj, test_ran_ef = test_ran_ef)
    # Check if results are valid
    if (any(LRT_results[4, ] < 0)) {
      # which_tests_weird <- colnames(LRT_results)[which(LRT_results[4, ] < 0)]
      which_tests_weird <- which(LRT_results[4, ] < 0)
      message(paste("Reduced model has a higher log likelihood for: ", paste(which_tests_weird, collapse = ", "), ". Results cannot be trusted! Try fitting the model with a different optimizer or reducing the random-effects structure.",
                    sep = ""))
    }
    if (tests == "bootstrap") {
      pb_objects <- lapply(reduced_fits, function(fit_tmp) {
        boot_tmp <- compute_parametric_bootstrap_test(fit_obj, fit_tmp, dv = dv, data = data,
                                                      distribution = distribution,
                                                      mm = mm, nsim = nsim, cl = cl, control = control, seed = seed)
        return(boot_tmp)
      })
      boot_table <- sapply(pb_objects, function(x) { return(data.frame(x$test)[2, ]) })

      # Return an mesdt_test object
      to_return <- new_mesdt_test(list(
        "fit_obj" = mesdt_fit,
        "reduced_fits" = reduced_fits,
        "LRT_results" = t(LRT_results),
        "type" = type,
        "pb_test_results" = t(boot_table),
        "pb_objects" = pb_objects,
        "seed" = seed
      ))

    } else {
      to_return <- new_mesdt_test(list(
        "fit_obj" = mesdt_fit,
        "reduced_fits" = reduced_fits,
        "LRT_results" = t(LRT_results),
        "type" = type
      ))
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
          boot_tmp <- compute_parametric_bootstrap_test(fit_full, fit_tmp, data = data, distribution = distribution,
                                                        mm = mm, dv = dv, nsim = nsim, cl = cl, control = control,
                                                        seed = seed)
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
          boot_tmp <- compute_parametric_bootstrap_test(fit_full, fit_tmp, data = data, distribution = distribution,
                                                        mm = mm, dv = dv, nsim = nsim, cl = cl, control = control,
                                                        seed = seed)
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
    # Check if results are valid
    if (any(LRT_results[4, ] < 0)) {
      which_tests_weird <- colnames(LRT_results)[which(LRT_results[4, ] < 0)]
      message(paste("Reduced model has a higher log likelihood for: ", paste(which_tests_weird, collapse = ", "), ". Results cannot be trusted! Try fitting the full model with a different optimizer or reducing the random-effects structure.",
                    sep = ""))
    }

    reduced_fits <- c(reduced_fits_lambda, reduced_fits_mu)
    full_fits <- c(full_fits_lambda, full_fits_mu)

    if (tests == "bootstrap") {

      pb_objects <- c(pb_objects_lambda, pb_objects_mu)
      boot_table <- cbind(boot_table_lambda, boot_table_mu)
      # Return a new mesdt_test object
      to_return <- new_mesdt_test(list(
        "fit_obj" = mesdt_fit,
        "reduced_fits" = reduced_fits,
        "LRT_results" = t(LRT_results),
        "type" = type,
        "pb_test_results" = t(boot_table),
        "pb_objects" = pb_objects,
        "seed" = seed
      ))

    } else {
      # Return a new mesdt_test object
      to_return <- new_mesdt_test(list(
        "fit_obj" = mesdt_fit,
        "reduced_fits" = reduced_fits,
        "LRT_results" = t(LRT_results),
        "type" = type
      ))
    }
  }
  return(to_return)
}


fit_submodels <- function(formula_mu, formula_lambda, dv, data, mm, type = 3, distribution,
                          test_intercepts = F, rem_ran_ef = F, correlate_sdt_params = T,
                          cl = NULL, tests_discriminability = "all", tests_response_bias = "all",
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
      # print(range_lambda)
      # print(range_mu)
    } else {
      if (length(attr(mm[["lambda"]], "assign")) == 1) range_lambda <- c()
      else range_lambda <- list(1:max(attr(mm[["lambda"]], "assign")))
      if (length(attr(mm[["mu"]], "assign")) == 1) range_mu <- c()
      else range_mu <- list(1:max(attr(mm[["mu"]], "assign")))
    }

    # Update range of to-be-tested parameters if only some parameters are to be tested
    if (is.null(tests_discriminability)) {
      labels_mu <- NULL; test_mu_labels <- NULL; range_mu <- NULL;
    } else if (tests_discriminability != "all") {
      labels_mu <- attr(stats::terms.formula(lme4::nobars(formula_mu)), "term.labels")
      test_mu_labels <- attr(stats::terms.formula(tests_discriminability), "term.labels")
      if (! all_terms_in(test_mu_labels, labels_mu)) stop("Only parameters that are in the model can be tested.")
      range_mu <- ifelse(test_intercepts,
                         list(range_mu[[1]][c(1, which_terms_in(reference_terms = labels_mu,
                                                                terms_to_check = test_mu_labels) + 1)]),
                         list(range_mu[[1]][which_terms_in(reference_terms = labels_mu,
                                                           terms_to_check = test_mu_labels)]))
    }
    if (is.null(tests_response_bias)) {
      labels_lambda <- NULL; test_lambda_labels <- NULL; range_lambda <- NULL;
    } else if (tests_response_bias != "all") {
      labels_lambda <- attr(stats::terms.formula(lme4::nobars(formula_lambda)), "term.labels")
      test_lambda_labels <- attr(stats::terms.formula(tests_response_bias), "term.labels")
      if (! all_terms_in(test_lambda_labels, labels_lambda)) stop("Only parameters that are in the model can be tested.")
      range_lambda <- ifelse(test_intercepts,
                             list(range_lambda[[1]][c(1, which_terms_in(reference_terms = labels_lambda,
                                                                        terms_to_check = test_lambda_labels) + 1)]),
                             list(range_lambda[[1]][which_terms_in(reference_terms = labels_lambda,
                                                                   terms_to_check = test_lambda_labels)]))
    }
    #print(range_lambda)
    # print(range_mu)

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

    if (is.null(tests_response_bias)) full_formulas_lambda <- NULL
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

    if (is.null(tests_discriminability)) full_formulas_mu <- NULL
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
    names_lambda <- c("(Intercept)_lambda")
    names_mu <- c("(Intercept)_mu")
  }
  if (! rem_ran_ef) {
    if (length(attr(terms(lme4::nobars(formula_lambda)), "term.labels")) > 0) {
      if (! is.null(tests_response_bias)) {
        if (tests_response_bias != "all") names_lambda <- paste(c(names_lambda, test_lambda_labels), "lambda", sep = "_")
        else names_lambda <- c(names_lambda, paste(attr(terms(lme4::nobars(formula_lambda)), "term.labels"), "lambda", sep = "_"))
      }
    }

    if (! is.null(tests_discriminability)) {
      if (length(attr(terms(lme4::nobars(formula_mu)), "term.labels")) > 0) {
        if (tests_discriminability != "all") names_mu <- paste(c(names_mu, test_mu_labels), "mu", sep = "_")
        else names_mu <- c(names_mu, paste(attr(terms(lme4::nobars(formula_mu)), "term.labels"), "mu", sep = "_"))
      }
    }
    if (length(reduced_formulas_lambda) != 0) names(reduced_formulas_lambda) <- names_lambda
    if (length(reduced_formulas_mu) != 0) names(reduced_formulas_mu) <- names_mu
  }

  if (type == 2) {
    if (is.null(cl)) {
      # print("No cluster")

      full_fits <- lapply(c(full_formulas_lambda, full_formulas_mu),
                          fit_glmm, data = data, mm = mm, distribution = distribution, dv = dv,
                          control = control)
      # for reduced fits: check if the fit is already in the full_fits

      full_formulas_char <- sapply(c(full_formulas_lambda, full_formulas_mu), function(x) {as.character(x)[3]})

      reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), function(formula_tmp) {
        if (as.character(formula_tmp)[3] %in% full_formulas_char) {
          #print(paste("Copying ", formula_tmp, sep = ""))
          which_model_equal <- which(full_formulas_char == as.character(formula_tmp)[3])
          #print(which_model_equal)
          return(full_fits[[which_model_equal]])
        } else return(fit_glmm(formula_tmp, data, mm, distribution, dv, control))
      })
      all_fits <- c(full_fits, reduced_fits)
    } else {
      print("cluster")
      throwaway <- parallel::clusterExport(cl = cl,
                                           varlist = c("data", "mm", "fit_glmm"),
                                           env = environment())
      all_fits <- parallel::clusterApplyLB(cl,
                                           c(full_formulas_lambda,
                                             full_formulas_mu,
                                             reduced_formulas_lambda,
                                             reduced_formulas_mu),
                                           fit_glmm,
                                           data = data, mm = mm, distribution = distribution, dv = dv,
                                           control = control)
      if (options("mesdt.backend") == "lme4") all_fits <- unlist(all_fits)
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
      #print("No cluster")
      reduced_fits <- lapply(c(reduced_formulas_lambda, reduced_formulas_mu), fit_glmm,
                             data = data, mm = mm, distribution = distribution, dv = dv,
                             control = control)
      #print("reduced_fits done")
    } else {
      #print("cluster")
      throwaway <- parallel::clusterExport(cl = cl,
                                           varlist = c("data", "mm", "fit_glmm", "control"),
                                           env = environment())
      reduced_fits <- parallel::clusterApplyLB(cl,
                                               c(reduced_formulas_lambda, reduced_formulas_mu),
                                               fit_glmm,
                                               data = data, mm = mm, distribution = distribution,
                                               control = control)
      if (options("mesdt.backend") == "lme4") reduced_fits <- unlist(reduced_fits)
      names(reduced_fits) <- names(c(reduced_formulas_lambda, reduced_formulas_mu))
    }
    print(names(reduced_fits))
    return(reduced_fits)
  }
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
                                              distribution,
                                              nsim = 1000, cl = NULL, control = NULL,
                                              seed = NULL) {
  #if (is.null(seed)) seed <- sample(1:1e6, 1)

  do_pb <- function(x) {
    dat_tmp <- data
    # sim_dat_tmp <- stats::simulate(small_model)[[1]]
    sim_dat_tmp <- stats::simulate(small_model, seed = seed + x)[[1]]
    # Do refitting manually
    dat_tmp[[dv]] <- sim_dat_tmp
    sim_fit_full <- fit_glmm(stats::formula(large_model), dat_tmp, mm, distribution, dv, control)
    sim_fit_red <- fit_glmm(stats::formula(small_model), dat_tmp, mm, distribution, dv, control)
    return(-2 * (stats::logLik(sim_fit_red) - stats::logLik(sim_fit_full)))
  }

  if (is.null(cl)) {
    parallel <- F
    LRs_boot <- sapply(1:nsim, do_pb)
  } else {
    parallel <- T
    # load variables on cluster
    throwaway <- parallel::clusterExport(cl = cl,
                            varlist = c("data", "mm", "fit_glmm", "dv", "small_model", "large_model", "control", "seed"),
                            #varlist = c("data", "mm", "fit_glmm", "dv", "small_model", "large_model", "control"),
                            env = environment())
    if (options("mesdt.backend") == "glmmTMB") throwaway <- parallel::clusterCall(cl = cl, "require", package = "glmmTMB", character.only = T)
    # set seed on cluster
    LRs_boot <- parallel::clusterApplyLB(cl, 1:nsim, do_pb)
    if (options("mesdt.backend") == "lme4") LRs_boot <- unlist(LRs_boot)
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
