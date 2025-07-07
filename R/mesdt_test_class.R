# Constructor
new_mesdt_test <- function(input_list) {
  structure(input_list, class = "mesdt_test")
}

#' print method
#' TODO
#' @export
print.mesdt_test <- function(obj) {
  print.summary.mesdt_test(obj)
}


#' @export
summary.mesdt_test <- function(obj) {
  to_return <- list()
  to_return[["LRT_results"]] <- obj$LRT_results
  to_return[["type"]] <- obj$type
  if (! is.null(obj$pb_test_results)) {
    to_return[["pb_test_results"]] <- obj$pb_test_results
  }
  return(structure(to_return, class = "summary.mesdt_test"))
}


#' @export
print.summary.mesdt_test <- function(x) {
  # if bootstrap exists:
  if (! is.null(x$pb_test_results)) {
    meth_str <- paste("Type ", ifelse(x$type == 2, "II", "III"),
                      " parametric bootstrap tests \n\n", sep = "")
    cat(meth_str)
    cat("Discriminability: \n")
    which_mu <- grep("mu", rownames(x$pb_test_results))
    tests_d <- data.frame(x$pb_test_results[which_mu, ])
    rownames(tests_d) <- gsub("_mu", "", rownames(x$pb_test_results)[which_mu])
    tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
    print(tests_d)
    cat("\n")
    cat("Response Bias: \n")
    which_c <- grep("lambda", rownames(x$pb_test_results))
    tests_c <- data.frame(x$pb_test_results[which_c, ])
    rownames(tests_c) <- gsub("_lambda", "", rownames(x$pb_test_results)[which_c])
    tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
    print(tests_c)
  } else {
    meth_str <- paste("Type ", ifelse(x$type == 2, "II", "III"),
                      " likelihood ratio tests \n\n", sep = "")
    cat(meth_str)
    cat("Discriminability: \n")
    which_mu <- grep("mu", rownames(x$LRT_results))
    tests_d <- data.frame(x$LRT_results[which_mu, ])
    rownames(tests_d) <- gsub("_mu", "", rownames(x$LRT_results)[which_mu])
    tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
    tests_d[, 2] <- sapply(tests_d[, 2], round, 2)
    tests_d[, 4] <- sapply(tests_d[, 4], round, 2)
    print(tests_d)
    cat("\n")
    cat("Response Bias: \n")
    which_c <- grep("lambda", rownames(x$LRT_results))
    tests_c <- data.frame(x$LRT_results[which_c, ])
    rownames(tests_c) <- gsub("_lambda", "", rownames(x$LRT_results)[which_c])
    tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
    tests_c[, 2] <- sapply(tests_c[, 2], round, 2)
    tests_c[, 4] <- sapply(tests_c[, 4], round, 2)
    print(tests_c)
  }

  invisible()
}
