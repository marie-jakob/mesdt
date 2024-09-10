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
                                    mm = NULL, to_remove = NULL, remove_correlations = F) {

  # check if all fixed effects are removed
  if (!is.null(to_remove[["lambda"]]) & !is.null(to_remove[["mu"]])) {
    if (to_remove[["lambda"]] == Inf & to_remove[["mu"]] == Inf) {
      message("Model must contain fixed effects.")
      return()
    }
  }


  # get random grouping factors for mu and lambda
  rdm_facs_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })
  rdm_facs_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
  })

  # Handle given random effects-structure and convert to the corresponding internal formula
  rdm_pred_lambda <- sapply(lme4::findbars(formula_lambda), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })
  rdm_pred_mu <- sapply(lme4::findbars(formula_mu), function(x) {
    gsub(" ", "", gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2]))
  })

  # add "_rdm_fac" to the model matrix to allow for multiple random effects
  if (! is.null(to_remove)) {
    tbr_char <- c()
    for (i in to_remove) {
      tmp <- deparse(i)
      if (grepl(":", tmp)) tmp <- paste("c(", tmp, ")", sep = "")
      tbr_char <- c(tbr_char, tmp)
    }
    names(tbr_char) <- names(to_remove)
  }

  # Case 1: Everything Correlated
  # -> one random-effects term
  # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)
  # 0 to suppress automatic intercept (is contained in the mm)
  rdm_formula_parts <- list()
  # Correlated random effects
  if (! remove_correlations & length(rdm_facs_lambda) == length(unique(rdm_facs_lambda)) &
      length(rdm_facs_mu) == length(unique(rdm_facs_mu))) {
    print("correlate_sdt_params")
    for (rdm_fac in rdm_facs_lambda) {
      name_lambda_tmp <- paste("rdm_lambda_", rdm_fac, sep = "")
      name_mu_tmp <- paste("rdm_mu_", rdm_fac, sep = "")
      # Three cases:
      # 1. add the whole matrix
      lambda_tmp <- ifelse(is.null(to_remove[[name_lambda_tmp]]), paste("mm[['", name_lambda_tmp, "']]", sep = ""),
                           # 2. remove the whole matrix
                           ifelse(all(to_remove[[name_lambda_tmp]] == Inf), "",
                                  # 3. remove the given indices
                                  paste("mm[['", name_lambda_tmp, "']][, -", tbr_char[[name_lambda_tmp]], "]", sep = "")
                                  ))
      if (rdm_fac %in% rdm_facs_mu) {
        # Three cases:
        # 1. add the whole matrix
        mu_tmp <- ifelse(is.null(to_remove[[name_mu_tmp]]), paste("mm[['", name_mu_tmp, "']]", sep = ""),
                             # 2. remove the whole matrix
                             ifelse(all(to_remove[[name_mu_tmp]] == Inf), "",
                                    # 3. remove the given indices
                                    paste("mm[['", name_mu_tmp, "']][, -", tbr_char[[name_mu_tmp]], "]", sep = "")
                             ))
        } else mu_tmp <- ""
      if (lambda_tmp == "" & mu_tmp == "") {
        form_part_tmp <- ""
      } else if (lambda_tmp == "") {
        form_part_tmp <- paste("(0 + ", mu_tmp, "| ", rdm_fac, ")", sep = "")
      } else if (mu_tmp == "") {
        form_part_tmp <- paste("(0 + ", lambda_tmp, "| ", rdm_fac, ")", sep = "")
      } else {
        if (correlate_sdt_params) form_part_tmp <- paste("(0 + ", paste(lambda_tmp, mu_tmp, sep = " + "), " | ", rdm_fac, ")", sep = "")
        else form_part_tmp <- paste("(0 + ", lambda_tmp, "| ", rdm_fac, ") + ", "(0 + ", mu_tmp, " | ", rdm_fac, ")", sep = "")
      }
      form_part_tmp[sapply(form_part_tmp, is.null)] <- NULL
      rdm_formula_parts <- append(rdm_formula_parts, form_part_tmp)
    }
    for (rdm_fac in rdm_facs_mu) {
      if (! rdm_fac %in% rdm_facs_lambda) {
        name_mu_tmp <- paste("rdm_mu_", rdm_fac, sep = "")
        # Three cases:
        # 1. add the whole matrix
        mu_tmp <- ifelse(is.null(to_remove[[name_mu_tmp]]), paste("mm[['", name_mu_tmp, "']]", sep = ""),
                          # 2. remove the whole matrix
                          ifelse(all(to_remove[[name_mu_tmp]] == Inf), "",
                                # 3. remove the given indices
                                paste("mm[['", name_mu_tmp, "']][, -", tbr_char[[name_mu_tmp]], "]", sep = "")
                           ))

        print(mu_tmp)
        if (mu_tmp != "") rdm_formula_parts <- append(rdm_formula_parts,
                                                           paste("(0 + ", mu_tmp, " | ", rdm_fac, ")", sep = ""))
      }
    }

    ################
    # Completely uncorrelated random effects
    ################
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
      name_lambda_tmp <- paste("rdm_lambda_", rdm_fac, sep = "")
      rdm_parts_tmp <- sapply(1:ncol(mm[[name_lambda_tmp]]), function(x) {
        if (x %in% to_remove[[name_lambda_tmp]]) return()
        else if (all(to_remove[[name_lambda_tmp]] == Inf) & ! is.null(to_remove[[name_lambda_tmp]])) return()
        else return(paste("(0 + mm[['rdm_lambda_", rdm_fac, "']][, ", x, "] | ", rdm_fac, ")", sep ="", recycle0 = T))
      })
      rdm_parts_tmp[sapply(rdm_parts_tmp, is.null)] <- NULL
      if (any(rdm_parts_tmp != "")) rdm_formula_parts <- append(rdm_formula_parts, paste(rdm_parts_tmp, collapse = "+"))

    }
    for (rdm_fac in unique(rdm_facs_mu)) {
      name_mu_tmp <- paste("rdm_mu_", rdm_fac, sep = "")
      rdm_parts_tmp <- sapply(1:ncol(mm[[name_mu_tmp]]), function(x) {
        if (x %in% to_remove[[name_mu_tmp]]) return()
        else if (all(to_remove[[name_mu_tmp]] == Inf) & ! is.null(to_remove[[name_mu_tmp]])) return()
        else return(paste("(0 + mm[['rdm_mu_", rdm_fac, "']][, ", x, "] | ", rdm_fac, ")", sep =""))
      })
      rdm_parts_tmp[sapply(rdm_parts_tmp, is.null)] <- NULL
      if (any(rdm_parts_tmp != "")) rdm_formula_parts <- append(rdm_formula_parts, paste(rdm_parts_tmp, collapse = "+"))
    }
  }
  rdm_formula <- paste(rdm_formula_parts, collapse = " + ")
  # Case 3: Mix -> TODO for later
  # -> some terms correlated, some not
  # -> ~ ... + (mm[["rdm_mu"]] + mm[["rdm_lambda]] | ID)

  # make the full model formula if no parameter index to be removed is given
  #fixed_formula <- "0 + " #mm[['lambda']] + mm[['mu']]"
  fixed_formula <- "0"

  if (! is.null(to_remove[["lambda"]])) {
    if (! all(to_remove[["lambda"]] == Inf)) fixed_formula <- paste(fixed_formula, paste("mm[['lambda']][, -", tbr_char[["lambda"]], ']', sep = ""), sep = " + ")
  } else fixed_formula <- paste(fixed_formula, "mm[['lambda']]", sep = " + ")

  if (! is.null(to_remove[["mu"]])) {
    if (! all(to_remove[["mu"]] == Inf)) fixed_formula <- paste(fixed_formula, paste("mm[['mu']][, -", tbr_char[["mu"]], ']', sep = ""), sep = " + ")
  } else fixed_formula <- paste(fixed_formula, "mm[['mu']]", sep = " + ")

  if ((is.null(lme4::findbars(formula_lambda)) & is.null(lme4::findbars(formula_mu))) |
      rdm_formula == "") {
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
