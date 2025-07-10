

# Constructor
new_mesdt_fit <- function(input_list) {
  structure(input_list, class = "mesdt_fit")
}


#' @export
print.mesdt_fit <- function(x) {
  print.summary.mesdt_fit(x)
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
    rownames(c_coef) <- substr(rownames(c_coef), 15, nchar(rownames(c_coef)))
  }
  if (obj$user_input$distribution != "gumbel-min") {
    c_coef[, 1:3] <- (-1) * c_coef[, 1:3]
  }

  # Prepare random effects
  if (obj$user_input$backend == "lme4") {
    cov_mat <- VarCorr(obj$fit_obj)
    mixed <- T
  } else if (obj$user_input$backend == "glmmTMB") {
    cov_mat <- VarCorr(obj$fit_obj)$cond
    mixed <- T
  } else mixed <- F

  if (mixed) {
    # if random effects are correlated, print one correlation matrix
    # Remove "mm" stuff from colnames

    for (fac in names(cov_mat)) {
      colnames_tmp <- colnames(cov_mat[[fac]])
      if (any(grepl("\\[, *[0-9]+\\]", colnames_tmp))) {
        # -> uncorrelated random effect: generate new name from model matrix column
        # only works for completely uncorrelated random effects like this
        col_num <- as.numeric(gsub(".*\\[, *([0-9]+)\\].*", "\\1", colnames_tmp))
        sdt_name <- ifelse(grepl("lambda", colnames_tmp), "lambda", "mu")
        ran_name <- paste("rdm", sdt_name, gsub("\\.*[0-9]+", "", fac), sep = "_")

        new_col_name <- paste(
          colnames(obj$internal$mm[[ran_name]])[col_num],
          ifelse(sdt_name == "lambda", "(Response Bias)", "(Discriminability)"), sep = "")
        colnames(cov_mat[[fac]]) <- new_col_name
        rownames(cov_mat[[fac]]) <- new_col_name

      } else {
        sdt_names <- ifelse(grepl("lambda", colnames_tmp), "(Response Bias)", "(Discriminability)")
        colnames_new <- colnames_tmp
        for (i in 1:length(colnames_new)) {
          colnames_new[i] <- paste(gsub('mm\\[\\[".*?"\\]\\]', '', colnames_tmp[i]),
                                   sdt_names[i], sep = "")
          if (colnames_new[i] %in% c("(Response Bias)", "(Discriminability)")) {
            mat_name <- sub('mm\\[\\["(.*?)"\\]\\].*', '\\1', colnames_tmp[i])
            colnames_new[i] <- paste(colnames(obj$internal$mm[[mat_name]]),
                                     sdt_names[i], sep = "")
          }
        }
        colnames(cov_mat[[fac]]) <- colnames_new
        rownames(cov_mat[[fac]]) <- colnames_new
      }
    }

    #fitMsgs <- lme4::.merMod.msgs(obj$fit_obj)
    #if(any(nchar(fitMsgs) > 0)) {
    #  cat("fit warnings:\n"); writeLines(fitMsgs)
    #}
    #.prt.warn(x@optinfo,summary=TRUE)

  }
  to_return <- list(
    "d_coef" = d_coef,
    "c_coef" = c_coef,
    "opt_info" = opt_info,
    "user_input" = obj$user_input
  )

  if (mixed) to_return[["cov_mat"]] <- cov_mat
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
    meth <- ifelse(nAGQ != 1, paste("(Adaptive Gauss-Hermite Quadrature, nAGQ = ", nAGQ, ")", sep = ""),
                   "(Laplace Approximation)")
    pr <- paste(pr, meth, sep = "")
  }
  pack <- x$user_input$backend
  if (pack != "glm") {
    pr <- paste(pr, " with the ", pack, " package. \n \n", sep = "")
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
                                    signif.stars = FALSE,
                                    ...) {

  printmethod(x)

  cat("Discriminability: ", deparse(x$user_input$discriminability), "\n")
  cat("Response Bias:    ", deparse(x$user_input$bias), "\n\n")

  # Print random effects
  if (! is.null(x$cov_mat)) {
    .prt.VC(x$cov_mat, digits)
    cat("\n")
  }


  # Print fixed effects
  cat("Fixed effects and Wald tests for discriminability: \n")
  stats::printCoefmat(x$d_coef, digits = digits, signif.stars = signif.stars)
  cat("\nFixed effects and Wald tests for response bias: \n")
  stats::printCoefmat(x$c_coef, digits = digits, signif.stars = signif.stars)

  if(length(x$fitMsgs) && any(nchar(x$fitMsgs) > 0)) {
    cat("fit warnings:\n"); writeLines(x$fitMsgs)
  }
  #.prt.warn(x$optinfo,summary=FALSE)
  #invisible(x)
}



#' @importFrom stats simulate
#' @export
simulate.mesdt_fit <- function(obj, ...) {
  # get method for correct backend
  #print(mesdt_obj$backend)
  #if (mesdt_obj$backend == "lme4") pred <- lme4::simulate.merMod(mesdt_obj$fit_obj, ...)
  #else if (mesdt_obj$backend == "glmmTMB") pred <- glmmTMB::simulate(mesdt_obj$fit_obj, ...)
  pred <- stats::simulate(obj$fit_obj, ...)
  return(pred)
}

#------------------------------------------------------------------------------#
#### From lme4 ####
# -> All of the following functions are taken from the lme4 package


##' "format()" the 'VarCorr' matrix of the random effects -- for
##' print()ing and show()ing
##'
##' @title Format the 'VarCorr' Matrix of Random Effects
##' @param varc a \code{\link{VarCorr}} (-like) matrix with attributes.
##' @param digits the number of significant digits.
##' @param comp character vector of length one or two indicating which
##' columns out of "Variance" and "Std.Dev." should be shown in the
##' formatted output.
##' @param formatter the \code{\link{function}} to be used for
##' formatting the standard deviations and or variances (but
##' \emph{not} the correlations which (currently) are always formatted
##' as "0.nnn"
##' @param ... optional arguments for \code{formatter(*)} in addition
##' to the first (numeric vector) and \code{digits}.
##' @return a character matrix of formatted VarCorr entries from \code{varc}.
formatVC <- function(varcor, digits = max(3, getOption("digits") - 2),
                     comp = "Std.Dev.", corr = any(comp == "Std.Dev."),
                     formatter = format,
                     useScale = attr(varcor, "useSc"),
                     ...) {
  c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
  avail.c <- c.nms[-(1:2)]
  if(anyNA(mcc <- pmatch(comp, avail.c)))
    stop("Illegal 'comp': ", comp[is.na(mcc)])
  nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
  if(length(use.c) == 0)
    stop("Must show variances and/or standard deviations")
  reStdDev <- c(lapply(varcor, attr, "stddev"),
                if(useScale) list(Residual = unname(attr(varcor, "sc"))))
  reLens <- lengths(reStdDev)
  nr <- sum(reLens)
  reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
  reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)
  reMat[,"Name"] <- c(unlist(lapply(varcor, colnames)), if(useScale) "")
  if (any("Variance" == use.c))
    reMat[,"Variance"] <- formatter(unlist(reStdDev)^2, digits = digits, ...)
  if (any("Std.Dev." == use.c))
    reMat[,"Std.Dev."] <- formatter(unlist(reStdDev),   digits = digits, ...)
  if (any(reLens > 1L)) { ## append lower triangular matrix of correlations / covariances
    maxlen <- max(reLens)
    Lcomat <- if(corr)
      lapply(varcor, attr, "correlation")
    else # just the matrix, i.e. {dim,dimnames}
      lapply(varcor, identity)
    ## function(v) `attributes<-`(v, attributes(v)[c("dim","dimnames")])
    co <- # corr or cov
      do.call(rbind,
              lapply(Lcomat,
                     function(x) {
                       x <- as.matrix(x)
                       dig <- max(2, digits - 2) # use 'digits' !
                       ## not using formatter() for correlations
                       cc <- format(round(x, dig), nsmall = dig)
                       cc[!lower.tri(cc)] <- ""
                       nr <- nrow(cc)
                       if (nr >= maxlen) return(cc)
                       cbind(cc, matrix("", nr, maxlen-nr))
                     }))[, -maxlen, drop = FALSE]
    if (nrow(co) < nrow(reMat))
      co <- rbind(co, matrix("", nrow(reMat) - nrow(co), ncol(co)))
    colnames(co) <- c(if(corr) "Corr" else "Cov", rep.int("", max(0L, ncol(co)-1L)))
    cbind(reMat, co, deparse.level=0L)
  } else reMat
}


.prt.VC <- function(varcor, digits,
                    comp = "Std.Dev.",
                    corr = any(comp == "Std.Dev."),
                    formatter = format, ...) { # '...' *only* passed to print()
  cat("Random effects:\n")
  fVC <- formatVC(varcor, digits=digits, formatter=formatter, comp=comp, corr=corr)
  print(fVC, quote = FALSE, digits = digits, ...)
}


.prt.warn <- function(optinfo, summary=FALSE, ...) {
  if(length(optinfo) == 0) return() # currently, e.g., from refitML()
  ## check all warning slots: print numbers of warnings (if any)
  cc <- optinfo$conv$opt
  msgs <- unlist(optinfo$conv$lme4$messages)
  ## can't put nmsgs/nwarnings compactly into || expression
  ##   because of short-circuiting
  nmsgs <- length(msgs)
  warnings <- optinfo$warnings
  nwarnings <- length(warnings)
  if (cc > 0 || nmsgs > 0 || nwarnings > 0) {
    m <- if (cc==0) {
      "(OK)"
    } else if (!is.null(optinfo$message)) {
      sprintf("(%s)",optinfo$message)
    } else ""
    convmsg <- sprintf("optimizer (%s) convergence code: %d %s",
                       optinfo$optimizer, cc, m)
    if (summary) {
      cat(convmsg,sprintf("; %d optimizer warnings; %d lme4 warnings",
                          nwarnings,nmsgs),"\n")
    } else {
      cat(convmsg,
          msgs,
          unlist(warnings),
          sep="\n")
      cat("\n")
    }
  }
}
