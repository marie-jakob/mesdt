# Constructor
new_mesdt_test <- function(input_list) {
  structure(input_list, class = "mesdt_test")
}

#' Print Method for mesdt_test Objects
#'
#' Prints a \code{mesdt_test} object by calling
#' \code{print.summary.mesdt_test()}.
#' @param x An object of class \code{mesdt_test}, typically returned by
#'   model comparison functions such as \code{mesdt()}.
#' @param ... Additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @return Invisibly returns \code{NULL}. The function prints results to the console.
#'
#' @method print mesdt_test
#' @export
#'
print.mesdt_test <- function(x, ...) {
  print.summary.mesdt_test(x)
}



#' Summary Method for mesdt_test Objects
#'
#' Provides a structured summary of a \code{mesdt_test} object, including
#' likelihood ratio test (LRT) results and (if
#' available) parametric bootstrap test results.
#'
#' @param object An object of class \code{mesdt_test}, typically returned by
#'   \code{mesdt_test()} or a related model comparison function.
#' @param ... Additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @return
#' A list of class \code{summary.mesdt_test} containing:
#' \describe{
#'   \item{\code{LRT_results}}{A table or list containing likelihood ratio
#'   test results.}
#'   \item{\code{type}}{A character string indicating the test type.}
#'   \item{\code{pb_test_results}}{Parametric bootstrap test results, if
#'   available.}
#' }
#'
#' @method summary mesdt_test
#' @export
#'
summary.mesdt_test <- function(object, ...) {
  obj <- object
  to_return <- list()
  to_return[["LRT_results"]] <- obj$LRT_results
  to_return[["type"]] <- obj$type
  if (! is.null(obj$pb_test_results)) {
    to_return[["pb_test_results"]] <- obj$pb_test_results
  }
  return(structure(to_return, class = "summary.mesdt_test"))
}


#' Print Method for Summary of mesdt_test Objects
#'
#' Prints a structured summary of a \code{summary.mesdt_test} object to the console.
#' Depending on whether parametric bootstrap results are available, it prints
#' either the parametric bootstrap tests or likelihood ratio test (LRT) results.
#'
#' @param x An object of class \code{summary.mesdt_test}, typically returned
#'   by \code{summary.mesdt_test()}.
#' @param ... Additional arguments passed to or from other methods. Currently unused.
#'
#' @return Invisibly returns \code{NULL}. The function prints results to the console.
#'
#' @export
print.summary.mesdt_test <- function(x, ...) {
  # if bootstrap exists:
  if (! is.null(x$pb_test_results)) {
    meth_str <- paste("Type ", ifelse(x$type == 2, "II", "III"),
                      " parametric bootstrap tests \n\n", sep = "")
    cat(meth_str)

    which_mu <- grep("mu", rownames(x$pb_test_results))
    if (! length(which_mu) == 0) {
      cat("Discriminability: \n")
      tests_d <- data.frame(x$pb_test_results[which_mu, ])
      rownames(tests_d) <- gsub("_mu", "", rownames(x$pb_test_results)[which_mu])
      tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
      print(tests_d)
      cat("\n")
    }
    which_c <- grep("lambda", rownames(x$pb_test_results))
    if (! length(which_c) == 0) {
      tests_c <- data.frame(x$pb_test_results[which_c, ])
      rownames(tests_c) <- gsub("_lambda", "", rownames(x$pb_test_results)[which_c])
      tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
      print(tests_c)
    }
  } else {
    meth_str <- paste("Type ", ifelse(x$type == 2, "II", "III"),
                      " likelihood ratio tests \n\n", sep = "")
    cat(meth_str)

    which_mu <- grep("mu", rownames(x$LRT_results))
    if (! length(which_mu) == 0) {

      cat("Discriminability: \n")

      tests_d <- data.frame(x$LRT_results[which_mu, ])
      rownames(tests_d) <- gsub("_mu", "", rownames(x$LRT_results)[which_mu])
      tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
      tests_d[, 2] <- sapply(tests_d[, 2], round, 2)
      tests_d[, 4] <- sapply(tests_d[, 4], round, 2)
      print(tests_d)
      cat("\n")
    }
    which_c <- grep("lambda", rownames(x$LRT_results))
    if (! length(which_c) == 0) {

      cat("Response Bias: \n")
      tests_c <- data.frame(x$LRT_results[which_c, ])
      rownames(tests_c) <- gsub("_lambda", "", rownames(x$LRT_results)[which_c])
      tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
      tests_c[, 2] <- sapply(tests_c[, 2], round, 2)
      tests_c[, 4] <- sapply(tests_c[, 4], round, 2)
      print(tests_c)
    }

  }

  invisible()
}
