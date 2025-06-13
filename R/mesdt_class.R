

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


  if (obj$user_input$backend == "lme4" | obj$user_input$backend == "glm") {

    # Get coefficients
    summ_lme4 <- summary(obj$fit_obj)
    d_coef <- summ_lme4$coefficients[grepl("mu", rownames(summ_lme4$coefficients)), ]
    c_coef <- summ_lme4$coefficients[grepl("lambda", rownames(summ_lme4$coefficients)), ]

    # Get optimizer info
    if (obj$user_input$backend == "lme4") {
      opt_info <- list(
        "nAGQ" = obj$fit_obj@devcomp$dims[["nAGQ"]]
      )
    } else opt_info <- NULL

  } else if (obj$user_input$backend == "glmmTMB") {
    summ_glmmtmb <- summary(obj$fit_obj)

    d_coef <- summ_glmmtmb$coefficients$cond[grepl("mu", rownames(summ_glmmtmb$coefficients$cond)), ]
    c_coef <- summ_glmmtmb$coefficients$cond[grepl("lambda", rownames(summ_glmmtmb$coefficients$cond)), ]
    opt_info <- NULL
  }

  if (is.null(rownames(d_coef))) {
    col_nms <- names(d_coef)
    d_coef <- matrix(d_coef, nrow = 1)
    colnames(d_coef) <- col_nms
    rownames(d_coef) <- colnames(obj$internal$mm$mu)
  } else {
    rownames(d_coef) <- substr(rownames(d_coef), 11, nchar(rownames(d_coef)))
  }
  if (is.null(rownames(c_coef))) {
    col_nms <-  names(c_coef)
    c_coef <- matrix(c_coef, nrow = 1)
    colnames(c_coef) <- col_nms
    rownames(c_coef) <- colnames(obj$internal$mm$lambda)
  } else {
    c_coef[, 1:3] <- (-1) * c_coef[, 1:3]
    rownames(c_coef) <- substr(rownames(c_coef), 15, nchar(rownames(c_coef)))
  }


  #fitMsgs <- lme4::.merMod.msgs(obj$fit_obj)
  #if(any(nchar(fitMsgs) > 0)) {
  #  cat("fit warnings:\n"); writeLines(fitMsgs)
  #}
  #.prt.warn(x@optinfo,summary=TRUE)

  to_return <- list(
    "user_input" = obj$user_input,
    "d_coef" = d_coef,
    "c_coef" = c_coef,
    "opt_info" = opt_info
  )
  if (! is.null(obj$LRTs)) to_return[["LRTs"]] <- obj$LRTs$LRT_results
  else if (! is.null(obj$PB_tests)) to_return[["PB_tests"]] <- obj$LRTs$PB_test_results

  return(structure(to_return, class = "summary.mesdt_fit"))
}


# Taken from lme4
printmethod <- function(x) {
  distr_pretty <- ifelse(x$user_input$distribution == "gaussian", "Gaussian",
                         ifelse(x$user_input$distr == "logistic", "logistic", "Gumbel-Min"))
  if (x$user_input$backend != "glm") pr <- "Mixed-effects signal "
  else pr <- "Signal "
  pr <- paste(pr, "detection theory model with ", distr_pretty, " evidence distributions fit by maximum likelihood ", sep = "")
  if (! is.null(x$opt_info)) {
    nAGQ <- x$opt_info$nAGQ
    meth <- ifelse(nAGQ != 1, paste(" (Adaptive Gauss-Hermite Quadrature, nAGQ = ", nAGQ, ")", sep = ""),
                   "(Laplace Approximation) ")
    pr <- paste(pr, meth, sep = "")
  }
  pack <- x$user_input$backend
  if (pack != "glm") {
    pr <- paste(pr, "with the ", pack, " package. \n \n", sep = "")
  } else {
    pr <- paste(pr, "with glm(). \n\n", sep = "")
  }
  cat(pr)
  return()
}


# print summary method
#' @export
print.summary.mesdt_fit <- function(x,
                                    digits = max(3, getOption("digits") - 3),
                                    signif.stars = FALSE) {


  printmethod(x)

  cat("Discriminability:", deparse(x$user_input$discriminability), "\n")
  cat("Response Bias:     ", deparse(x$user_input$bias), "\n\n")

  # Print random effects
  # cor_mat <- VarCorr(x$fit_obj)
  #if (x$user_input$correlate_sdt_params = T) {
    # if random effects are correlated, print one correlation matrix

  #} else {
    # two separate matrices
  #}



  # Print fixed effects
  cat("Fixed effects and Wald tests for discriminability: \n")
  stats::printCoefmat(x$d_coef, digits = digits, signif.stars = signif.stars)
  cat("\nFixed effects and Wald tests for response bias: \n")
  stats::printCoefmat(x$c_coef, digits = digits, signif.stars = signif.stars)

  if(length(x$fitMsgs) && any(nchar(x$fitMsgs) > 0)) {
    cat("fit warnings:\n"); writeLines(x$fitMsgs)
  }
  .prt.warn(x$optinfo,summary=FALSE)
  invisible(x)
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


