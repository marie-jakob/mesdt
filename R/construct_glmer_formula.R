#' Convert two given formulas for the two SDT parameters to a glmer() formula
#'
#' @param formula_mu formula for sensitivity
#' @param formula_lambda formula for response bias
#' @param dv name of dependent variable
#' @param correlate_sdt_params model correlations between mu and lambda random effects?
#' @param mm model matrices (corresponding to the formulas) -> necessary for formulas with uncorrelated
#'    random effects and for indices to be removed
#' @param param_idc optional vector of parameters indices to be removed to construct a reduced formula
#' @param remove_from_mu optional argument to indicate whether the to-be-removed parameter
#'    should be removed from mu or lambda model matrix
#'
#' @return lme4 formula
#'
#' @importFrom stats formula
#' @importFrom lme4 findbars
#'
#' @examples
#' construct_glmer_formula(
#'   formula_mu = ~ x1 + (1 | ID)
#'   formula_lambda = ~ x2 + (1 | ID)
#'   dv = "y",
#'   trial_type_var = "trial_type"
#' )
construct_glmer_formula <- function(formula_mu, formula_lambda, dv, correlate_sdt_params = T,
                                    mm = NULL, param_idc = NULL, remove_from_mu = F,
                                    remove_from_rdm = "") {

  # check if the random grouping factor is the same for mu and lambda
  rdm_facs_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })
  rdm_facs_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })
  # Handle given random effects-structure and convert to the corresponding internal formula
  # Logic: Variables are in the same parentheses if they are correlated (no "||" stuff)

  rdm_pred_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })
  rdm_pred_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })
  # add "_rdm_fac" to the model matrix to allow for multiple random effects
  # TODO: adapt for multiple factors

  if (! is.null(param_idc)) {
    param_idc_char <- deparse(param_idc)
    if (grepl(":", param_idc_char)) {
      param_idc_char <- paste("c(", param_idc_char, ")", sep = "")
    }
  }

  # Case 1: Everything Correlated
  # -> one random-effects term
  # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)
  # 0 to suppress automatic intercept (is contained in the mm)
  rdm_formula_parts <- list()
  if (length(rdm_facs_lambda) == length(unique(rdm_facs_lambda)) &
      length(rdm_facs_mu) == length(unique(rdm_facs_mu))) {
    if (correlate_sdt_params) {
      for (rdm_fac in rdm_facs_lambda) {
        if (rdm_fac %in% rdm_facs_mu) {
          if (remove_from_rdm == rdm_fac) {
            if (remove_from_mu) {
              form_part_tmp <- ifelse(all(param_idc == Inf),
                                      # remove whole random mu model matrix
                                      paste("(0 + mm[['rdm_lambda_", rdm_fac, "']] | ", rdm_fac, ")", sep = ""),
                                      paste("(0 + mm[['rdm_lambda_", rdm_fac, "']] + mm[['rdm_mu_", rdm_fac, "']][, -", param_idc_char, "] | ", rdm_fac, ")", sep = ""))
            } else {
              form_part_tmp <- ifelse(all(param_idc == Inf),
                                      # remove whole random lambda model matrix
                                      paste("(0 + mm[['rdm_mu_", rdm_fac, "']] | ", rdm_fac, ")", sep = ""),
                                      paste("(0 + mm[['rdm_lambda_", rdm_fac, "']][, -", param_idc_char, "] + mm[['rdm_mu_", rdm_fac, "']] | ", rdm_fac, ")", sep = ""))
            }

          } else {
            form_part_tmp <- paste("(0 + mm[['rdm_lambda_", rdm_fac, "']] + mm[['rdm_mu_", rdm_fac, "']] | ", rdm_fac, ")", sep = "")
          }
        } else {
          form_part_tmp <- paste("(0 + mm[['rdm_lambda_", rdm_fac, "']] | ", rdm_fac, ")", sep = "")
        }
        rdm_formula_parts <- append(rdm_formula_parts, form_part_tmp)
      }
      for (rdm_fac in rdm_facs_mu) {
        if (! rdm_fac %in% rdm_facs_lambda) {
          if (remove_from_rdm == rdm_fac & remove_from_mu) {
            form_part_tmp <- paste("(0 + mm[['rdm_mu_", rdm_fac, "']][, -", param_idc_char, "] | ", rdm_fac, ")", sep = "")
          } else {
            form_part_tmp <- paste("(0 + mm[['rdm_mu_", rdm_fac, "']] | ", rdm_fac, ")", sep = "")
          }
          rdm_formula_parts <- append(rdm_formula_parts, form_part_tmp)
        }
      }
    } else {
      for (rdm_fac in rdm_facs_lambda) {
        if (remove_from_rdm == rdm_fac & ! remove_from_mu) {
          rdm_formula_parts <- append(rdm_formula_parts,
                                      paste("(0 + mm[['rdm_lambda_", rdm_fac, "']][, -", param_idc_char, "] | ", rdm_fac, ")", sep = ""))
        } else {
          rdm_formula_parts <- append(rdm_formula_parts,
                                      paste("(0 + mm[['rdm_lambda_", rdm_fac, "']] | ", rdm_fac, ")", sep = ""))
        }

      }
      for (rdm_fac in rdm_facs_mu) {
        if (remove_from_rdm == rdm_fac & remove_from_mu) {
          rdm_formula_parts <- append(rdm_formula_parts,
                                      paste("(0 + mm[['rdm_mu_", rdm_fac, "']][, -", param_idc_char, "] | ", rdm_fac, ")", sep = ""))
          } else {
            rdm_formula_parts <- append(rdm_formula_parts,
                                        paste("(0 + mm[['rdm_mu_", rdm_fac, "']] | ", rdm_fac, ")", sep = ""))
          }
      }
    }
  } else {
    message("Given random-effects structure contains uncorrelated terms. Modeling all random effects
            parametes as uncorrelated since a mix of correlated and uncorrelated terms is not
            supported at the moment.")
    if (correlate_sdt_params) {
      message("Correlating SDT Parameters is not possible in the presence of uncorrelated terms.")
    }
    if (is.null(mm)) {
      message("Model matrices are needed to generate a formula for uncorrelated random terms.")
      return()
    }
    # Case 2: No Correlations
    # -> as many random-effects terms as predictors in the model

    for (rdm_fac in unique(rdm_facs_lambda)) {
      rdm_formula_parts <- append(rdm_formula_parts, paste(sapply(1:ncol(mm[[paste("rdm_lambda_", rdm_fac, sep = "")]]), function(x) {
        if (rdm_fac == remove_from_rdm & ! remove_from_mu & x %in% param_idc) return()
        else return(paste("(0 + mm[['rdm_lambda_", rdm_fac, "']][, ", x, "] | ", rdm_fac, ")", sep ="", recycle0 = T))
      }), collapse = "+"))
    }
    for (rdm_fac in unique(rdm_facs_mu)) {
      rdm_formula_parts <- append(rdm_formula_parts, paste(sapply(1:ncol(mm[[paste("rdm_mu_", rdm_fac, sep = "")]]), function(x) {
        if (rdm_fac == remove_from_rdm & remove_from_mu & x %in% param_idc) return()
        else return(paste("(0 + mm[['rdm_mu_", rdm_fac, "']][, ", x, "] | ", rdm_fac, ")", sep =""))
      }), collapse = "+"))
    }

  }
  rdm_formula <- paste(rdm_formula_parts, collapse = " + ")
  # Case 3: Mix -> TODO for later
  # -> some terms correlated, some not
  # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)

  # make the full model formula if no parameter index to be removed is given
  fixed_formula <- "0 + mm[['lambda']] + mm[['mu']]"
  if (remove_from_rdm == ""&  !is.null(param_idc) ) {
    # check if everything is removed -> remove the whole matrix then
    if (!is.null(param_idc) & all(param_idc == Inf)) {
      if (remove_from_mu) fixed_formula <- "0 + mm[['lambda']]"
      else fixed_formula <- "0 + mm[['mu']]"
    } else {
      # make a reduced model formula
      if (remove_from_mu) {
        #print(deparse(param_idc))
        # deparse() vector for nonstandard evaluation
        term_reduced <- paste("mm[['mu']][, -", param_idc_char, ']', sep = "")
        fixed_formula <- paste("0 + mm[['lambda']] + ", term_reduced, sep = "")
      } else {
        term_reduced <- paste("mm[['lambda']][, -", param_idc_char, ']', sep = "")
        fixed_formula <- paste("0 + ", term_reduced, " + mm[['mu']]", sep = "")
      }
    }

  }
  if (is.null(lme4::findbars(formula_lambda)) & is.null(lme4::findbars(formula_mu))) {
    glmer_formula <- as.formula(paste(dv, " ~ ", fixed_formula, sep = ""),
                                # parent.frame() sets scope of the parent environment (i.e., where the function
                                # is called from) for the formula -> necessary such that model matrices can be found
                                env = parent.frame())

  } else {
    glmer_formula <- as.formula(paste(dv, " ~ ", fixed_formula, " + ", rdm_formula, sep = ""),
                                # parent.frame() sets scope of the parent environment (i.e., where the function
                                # is called from) for the formula -> necessary such that model matrices can be found
                                env = parent.frame())
  }

  return(glmer_formula)
}
