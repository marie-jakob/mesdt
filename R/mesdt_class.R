



# Constructor
new_mesdt_fit <- function(input_list) {
  structure(input_list, class = "mesdt_fit")
}



#' print method
#' TODO
#' @export
#print.mesdt_fit <- function(obj) {
#  print("this is an mesdt object")
#}
print.mesdt_fit <- function(obj) {
  cat("Printing mesdt_fit object:\n")
  print(obj$data)
}



#' summary method
#' TODO
#' @export
summary.mesdt_fit <- function(obj) {


  if (obj$user_input$backend == "lme4") {

    # Get coefficients
    summ_lme4 <- summary(obj$fit_obj)
    d_coef <- summ_lme4$coefficients[grepl("lambda", rownames(summ_lme4$coefficients)), ]
    c_coef <- summ_lme4$coefficients[grepl("mu", rownames(summ_lme4$coefficients)), ]

    # Get optimizer info
    opt_info <- list(
      "nAGQ" = obj$fit_obj@devcomp$dims[["nAGQ"]]
    )

  } else if (obj$user_input$backend == "glmmTMB") {
    summ_glmmtmb <- summary(obj$fit_obj)

    d_coef <- summ_glmmtmb$coefficients$cond[grepl("lambda", rownames(summ_glmmtmb$coefficients$cond)), ]
    c_coef <- summ_glmmtmb$coefficients$cond[grepl("mu", rownames(summ_glmmtmb$coefficients$cond)), ]
    # Get optimizer info
    opt_info <- NULL
  }


  to_return <- list(
    "user_input" = obj$user_input,
    "d_coef" = d_coef,
    "c_coef" = c_coef,
    "opt_info" = opt_info
  )
  if (! is.null(obj$LRTs)) to_return[["LRTs"]] <- obj$LRTs$LRT_results
  else if (! is.null(obj$PB_tests)) to_return[["PB_tests"]] <- obj$LRTs$PB_test_results

  return(structure(to_return, class = "summary.mesdt_fit"))

  #mm <- obj$mm

  # Post-Processing the lme4 output
  # backend = ifelse(options("backend") == "", "lme4", options("backend"))
  #if (obj$backend == "lme4") {
  #  coefs_lambda <- summary(obj$fit_obj)$coefficients[grepl("lambda", rownames(summary(obj$fit_obj)$coefficients)), ]
  #  coefs_mu <- summary(obj$fit_obj)$coefficients[grepl("mu", rownames(summary(obj$fit_obj)$coefficients)), ]
  #} else if (obj$backend == "glmmTMB") {
  #  coefs_lambda <- summary(obj$fit_obj)$coefficients$cond[grepl("lambda", rownames(summary(obj$fit_obj)$coefficients$cond)), ]
  #  coefs_mu <- summary(obj$fit_obj)$coefficients$cond[grepl("mu", rownames(summary(obj$fit_obj)$coefficients$cond)), ]
  #}
  #rownames(coefs_lambda) <- gsub('mm', "", rownames(coefs_lambda))
  #rownames(coefs_lambda) <- colnames(mm[["lambda"]])
  #if (is.null(nrow(coefs_lambda))) {
  #  coefs_lambda <- t(data.frame(coefs_lambda))
  #} else {
  #  coefs_lambda <- data.frame(coefs_lambda)
  #}
  #rownames(coefs_lambda) <- colnames(mm[["lambda"]])

  #if (is.null(nrow(coefs_mu))) {
  #  coefs_mu <- t(data.frame(coefs_mu))
  #} else {
  #  coefs_mu <- data.frame(coefs_mu)
  #}
  #rownames(coefs_mu) <- colnames(mm[["mu"]])

  #print(coefs_mu)
  #print(coefs_lambda)
  #return(list("mu" = coefs_mu,
  #            "lambda" = coefs_lambda))
}



printmethod <- function(x) {
  pr <- "Mixed-effects signal detection theory model fit by maximum likelihood"
  if (! is.null(x$opt_info)) {
    nAGQ <- x$opt_info$nAGQ
    print(nAGQ)
    meth <- ifelse(nAGQ != 1, paste(" (Adaptive Gauss-Hermite Quadrature, nAGQ = ", nAGQ, ")", sep = ""),
                   "(Laplace Approximation) ")
    pr <- paste(pr, meth, sep = "")
  }
  pack <- x$mesdt_fit$user_input$backend
  pr <- paste(pr, " with the ", pack, " package. \n \n", sep = "")
  cat(pr)
  return()
}


# print summary method
#' @export
print.summary.mesdt_fit <- function(x,
                                    digits = max(3, getOption("digits") - 3),
                                    signif.stars = FALSE,
                                    aic = FALSE) {


  printmethod(x)

  cat("Formula sensitivity:   ", deparse(x$user_input$formula_mu), "\n")
  cat("Formula response bias: ", deparse(x$user_input$formula_lambda), "\n\n")


  # TODO: print AIC and log likelihood if requested
  # if ()
  # Print random effects


  # Print fixed effects
  cat("Fixed effects and Wald tests for sensitivity: \n")
  stats::printCoefmat(x$d_coef, digits = digits, signif.stars = signif.stars)
  cat("\nFixed effects and Wald tests for response bias: \n")
  stats::printCoefmat(x$c_coef, digits = digits, signif.stars = signif.stars)


  # if the model object has LRTs, print those:
  if (! is.null(x$LRTs)) {
    cat("\nLikelihood ratio tests (type X) for sensitivity: \n \n")
    LRTs_round <- x$LRTs
    LRTs_round[] <- lapply(x$LRTs, round, 3)
    print(LRTs_round)
  }

  # if the model object has LRTs, print those:
  if (! is.null(x$PB_tests)) {
    cat("\nParametric bootstrapping tests (type X) for sensitivity: \n \n")
    LRTs_round <- x$LRTs
    LRTs_round[] <- lapply(x$LRTs, round, 3)
    print(LRTs_round)
  }

}



# TODO: test if this works
#' @export
simulate.mesdt_fit <- function(mesdt_obj, ...) {
  # get method for correct backend
  #print(mesdt_obj$backend)
  #if (mesdt_obj$backend == "lme4") pred <- lme4::simulate.merMod(mesdt_obj$fit_obj, ...)
  #else if (mesdt_obj$backend == "glmmTMB") pred <- glmmTMB::simulate(mesdt_obj$fit_obj, ...)
  pred <- stats::simulate(mesdt_obj$fit_obj, ...)
  return(pred)
}


