#' Convert two given formulas for the two SDT parameters to an lme4 formula
#'
#' @param form_mu formula for sensitivity
#' @param form_lambda formula for response bias
#' @param dv name of dependent variable
#' @param trial_type_var name of variable coding the type of trial (signal vs. noise)
#' @param random name of the variable of the random factor
#' @param between predictors varying between levels of the random factor
#' @param within predictors varying within levels of the random factor
#'
#' @return lme4 formula
#' @export
#'
#' @importFrom stats formula
#'
#' @examples
#' make_glmer_formula(
#'   form_mu = mu ~ x_test,
#'   form_lambda = lambda ~ x_test,
#'   dv = "Y",
#'   trial_type_var = "trial_type",
#'   random = "ID",
#'   within = c("x_test")
#'  )
make_glmer_formula <- function(form_mu, form_lambda, dv,
                             trial_type_var = "trial_type",
                             random,
                             between = NULL,
                             within = NULL) {

  # TODO: add catches for submitted variables

  # extract predictors for sensitivity and response bias from formulas
  pred_mu <- all.vars(form_mu)[[2:length(all.vars(form_lambda))]]
  pred_lambda <- all.vars(form_lambda)[2:length(all.vars(form_lambda))]

  # build a new formula according to lme4 syntax
  # fixed effects
  fixed_terms <- paste(dv, "~", trial_type_var, sep = " ")

  for (pred in pred_lambda) fixed_terms <- paste(fixed_terms, pred, sep = " + ")

  for (pred in pred_mu) {
    new_term <- paste(pred, ":", trial_type_var, sep = "")
    fixed_terms <- paste(fixed_terms, new_term, sep = " + ")
  }

  # random effects
  random_terms <- paste("(", trial_type_var, sep = "")

  for (i in within) {
    if (i %in% pred_lambda) random_terms <- paste(random_terms, i, sep = " + ")
    if (i %in% pred_mu) {
      new_term <- paste(pred, ":", trial_type_var, sep = "")
      random_terms <- paste(random_terms, new_term, sep = " + ")
    }
  }
  random_terms <- paste(random_terms, " | ", random, ")", sep = "")

  form_sdt <- as.formula(paste(fixed_terms, random_terms, sep = " + "),
                         env = globalenv())

  return(form_sdt)
}
