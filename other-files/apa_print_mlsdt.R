# Formats a given glmer() multi-level SDT object according to APA standards to
# include in an RMarkdown manuscript
# Args:
#
#
#
#     coef_pretty: Subscript of the coefficient for the formatted output
# Returns: Apa-formatted string
#' Fit a multilevel signal detection theory model
#'
#' @param fit fitted glmer() model
#' @param tests object returned by compute_tests(), including the parameter in question
#' @param coef internal name of the coefficient (name in the glmer fit)
#' @param coef_pretty Subscript of the coefficient for the formatted output
#'
#' @return list
#' @importFrom lme4 fixef
#' @importFrom stats qnorm
#' @importFrom lme4 VarCorr
#' @importFrom stats vcov
#' @export
#'
#' @examples
apa_print_mlsdt <- function(fit, tests, coef, coef_pretty = NULL,
                            direction = NULL, CI = T) {
  fit <- fit$fit_obj
  # Get actual numbers
  fixed_ef <- ifelse(grepl("status", coef), 2 * lme4::fixef(fit)[[coef]], (-1) * lme4::fixef(fit)[[coef]])
  model_table <- summary(fit)$coef
  SE <- model_table[which(names(model_table[, 1]) == coef), 2]
  if (!is.null(direction)) {
    if (direction == "greater") {
      lower <- fixed_ef - stats::qnorm(0.95) * SE
      CI_string <- paste("95\\% CI $[", printnum(lower), ", ", "\\infty", "]$, ", sep = "")
    } else if (direction == "smaller") {
      upper <- fixed_ef + stats::qnorm(0.95) * SE
      CI_string <- paste("95\\% CI $[", "\\infty", ", ", papaja::printnum(upper), "]$, ", sep = "")
    }
  } else {
    lower <- fixed_ef - stats::qnorm(0.975) * SE
    upper <- fixed_ef + stats::qnorm(0.975) * SE
    CI_string <- paste("95\\% CI $[", papaja::printnum(lower), ", ", papaja::printnum(upper), "]$, ", sep = "")
  }
  if (! is.null(coef_pretty)) beta <- paste("$\\beta_{", coef_pretty, "}^{", superscript, "} = ", printnum(fixed_ef), "$, ", sep = "")
  else beta <- paste("$\\hat{\\beta_{", superscript, "}} = ", papaja::printnum(fixed_ef), "$, ", sep = "")
  chi_sq <- paste("$\\chi^2(1) = ", papaja::printnum(LRTs$lrt_results$chi_sq[LRTs$lrt_results$param == coef]), "$, ", sep = "")
  p_val <- LRTs$lrt_results$p_value[LRTs$lrt_results$param == coef]
  if (! is.null(direction)) {
    p_val <- p_val / 2
    print(printnum(fixed_ef))
    if (direction == "smaller" & fixed_ef > 0) p_val <- 1 - p_val
    else if (direction == "greater" & fixed_ef < 0) p_val <- 1 - p_val
    p_letter <- "$p_{\\text{one-tailed}}"
    #if (direction == "smaller") p_letter <- "p_{< 0}"
    #else if (direction == "greater") p_letter <- "p_{> 0}"
  } else p_letter <- "$p"
  if (p_val < .001) p <- paste(p_letter, " ", printp(p_val), "$", sep = "")
  else if (p_val == 0.999) p <- paste(p_letter, " > ", p_val, "$", sep = "")
  else p <- paste(p_letter, " = ", printp(p_val), "$", sep = "")
  return(paste(beta, CI_string, chi_sq, p, sep = ""))
}
