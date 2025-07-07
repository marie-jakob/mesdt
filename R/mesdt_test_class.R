# Constructor
new_mesdt_test <- function(input_list) {
  structure(input_list, class = "mesdt_test")
}

#' print method
#' TODO
#' @export
print.mesdt_test <- function(obj) {
  # if bootstrap exists:
  if (! is.null(obj$pb_objects)) {
    meth_str <- paste("Type ", ifelse(obj$type == 2, "II", "III"),
                      "parametric bootstrap tests \n\n", sep = "")
    cat("Discriminability: \n")
    tests_d <- obj$LRT_results[grep("mu", rownames(obj$LRT_results)), ]
    rownames(tests_d) <- gsub("_mu", "", rownames(tests_d))
    tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
    tests_d[, 2] <- sapply(tests_d[, 2], round, 2)
    tests_d[, 4] <- sapply(tests_d[, 4], round, 2)
    print(tests_d)
    cat("\n")
    cat("Response Bias: \n")
    tests_c <- obj$LRT_results[grep("lambda", rownames(obj$LRT_results)), ]
    rownames(tests_c) <- gsub("_lambda", "", rownames(tests_c))
    tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
    tests_c[, 2] <- sapply(tests_c[, 2], round, 2)
    tests_c[, 4] <- sapply(tests_c[, 4], round, 2)
    print(tests_c)
  }

  meth_str <- paste("Type ", ifelse(obj$type == 2, "II", "III"),
                    "likelihood ratio tests \n\n", sep = "")
  cat(meth_str)
  cat("Discriminability: \n")
  tests_d <- obj$LRT_results[grep("mu", rownames(obj$LRT_results)), ]
  rownames(tests_d) <- gsub("_mu", "", rownames(tests_d))
  tests_d[, 1] <- sapply(tests_d[, 1], round, 2)
  tests_d[, 2] <- sapply(tests_d[, 2], round, 2)
  tests_d[, 4] <- sapply(tests_d[, 4], round, 2)
  print(tests_d)
  cat("\n")
  cat("Response Bias: \n")
  tests_c <- obj$LRT_results[grep("lambda", rownames(obj$LRT_results)), ]
  rownames(tests_c) <- gsub("_lambda", "", rownames(tests_c))
  tests_c[, 1] <- sapply(tests_c[, 1], round, 2)
  tests_c[, 2] <- sapply(tests_c[, 2], round, 2)
  tests_c[, 4] <- sapply(tests_c[, 4], round, 2)
  print(tests_c)
  invisible()
}

