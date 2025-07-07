.onLoad <- function(libname, pkgname) {
  # TODO: check if package loaded / installed
  message("Setting backend to lme4.")
  options(mesdt.backend = "lme4")
}

all_terms_in <- function(terms_to_check, reference_terms) {
  normalize_interaction <- function(term) {
    parts <- unlist(strsplit(term, ":"))
    paste(sort(parts), collapse = ":")  # Sort components and reassemble
  }

  # Normalize all terms
  normalized_to_check <- sapply(terms_to_check, normalize_interaction)
  normalized_reference <- sapply(reference_terms, normalize_interaction)

  # Find matches
  all(match(normalized_to_check, normalized_reference, nomatch = 0) > 0)
}


same_interactions <- function(f1, f2) {
  extract_interactions <- function(formula) {
    terms <- attr(terms(formula), "term.labels")
    interactions <- grep(":", terms, value = TRUE)
    interactions <- lapply(strsplit(interactions, ":"), sort)  # Sort within interactions
    interactions <- sapply(interactions, paste, collapse = ":")  # Convert back to string
    sort(interactions)  # Sort the entire list for comparison
  }

  i1 <- extract_interactions(f1)
  i2 <- extract_interactions(f2)

  identical(i1, i2)
}


which_terms_in <- function(reference_terms, terms_to_check) {
  normalize_interaction <- function(term) {
    parts <- unlist(strsplit(term, ":"))
    paste(sort(parts), collapse = ":")  # Sort components and reassemble
  }

  # Normalize all terms
  normalized_to_check <- sapply(terms_to_check, normalize_interaction)
  normalized_reference <- sapply(reference_terms, normalize_interaction)

  # Find matches
  matches <- match(normalized_to_check, normalized_reference, nomatch = 0)
  return(matches[matches > 0])
}




pchisqmix <- function(q, df, mix, lower.tail = TRUE) {
  df_vec <- rep(df, length(q))
  mix_vec <- rep(mix, length(q))
  upper <- stats::pchisq(q = q, df = df, lower.tail = lower.tail)
  lower <- ifelse(df_vec == 1, if (lower.tail) 1 else 0,
                  stats::pchisq(q, df-1, lower.tail = lower.tail))
  return(mix_vec * lower + (1 - mix_vec) * upper)
}




standardize_dist_input <- function(x) {
  if (x %in% c("Gaussian", "gaussian", "normal", "Normal")) dist_std <- "gaussian"
  else if (x %in% c("logistic", "Logistic")) dist_std <- "logistic"
  else if (x %in% c("gumbel-min", "gumbel_min", "Gumbel-Min", "Gumbel-min",
                    "Gumbel_Min", "Gumbel_min", "extreme-value-min", "extreme_value_min",
                    "Extreme-value-min", "Extreme_value_min")) dist_std <- "gumbel-min"
  else dist_std <- NULL
  return(dist_std)
}


standardize_type_input <- function(x) {
  if (x %in% c(2, "II", "ii", "2", "two")) type_std <- 2
  else if (x %in% c(3, "III", "iii", "3", "trhee")) type_std <- 3
  else type_std <- NULL
  return(type_std)
}


standardize_tests_input <- function(x) {
  if (x %in% c("LRT", "lrt")) test_std <- "lrt"
  else if (x %in% c("boot", "Boot", "bootstrap", "Bootstrap", "pb", "PB")) test_std <- "bootstrap"
  else test_std <- NULL
  return(test_std)
}



check_sensitivity <- function(fit_obj) {
  summ_mesdt <- summary(fit_obj)
  mu_mean <- summ_mesdt$d_coef[rownames(summ_mesdt$d_coef) == "(Intercept)", 1]
  if (length(mu_mean) > 0) {
    if (mu_mean < 0) warning("Mean Population Sensitivity is < 0, indicating that the trial_type variable might be coded reversely.")
  }
}


standardize_fit_formulas <- function(form_disc, form_bias) {
  # Check if formula has a left-hand side
  form_disc_std <- ""
  form_bias_std <- ""
  if (length(as.character(form_disc)) == 3) {
    if (as.character(form_disc)[2] == "discriminability") {
      form_disc_std <- as.formula(paste(as.character(form_disc)[c(1, 3)],
                                        collapse = " "),
                                  env = globalenv())
    } else if (as.character(form_disc)[2] == "bias") {
      if (as.character(form_bias)[2] == "discriminability") {
        form_bias_std <- as.formula(paste(as.character(form_disc)[c(1, 3)],
                                          collapse = " "),
                                    env = globalenv())
      } else {
        stop("Cannot interpret formula input. Please specify formulas as either
           discriminability = ~ x and bias = ~ x or discriminability ~ x and bias ~ x.")
      }
    } else {
      stop("Cannot interpret formula input. Please specify formulas as either
           discriminability = ~ x and bias = ~ x or discriminability ~ x and bias ~ x.")
    }
  } else form_disc_std <- form_disc

  if (length(as.character(form_bias)) == 3) {
    if (as.character(form_bias)[2] == "bias") {
      form_bias_std <- as.formula(paste(as.character(form_bias)[c(1, 3)],
                                        collapse = " "),
                                  env = globalenv())
    } else if (as.character(form_bias)[2] == "discriminability") {
      if (as.character(form_disc)[2] == "bias") {
        form_disc_std <- as.formula(paste(as.character(form_bias)[c(1, 3)],
                                          collapse = " "),
                                    env = globalenv())
      } else {
        stop("Cannot interpret formula input. Please specify formulas as either
           discriminability = ~ x and bias = ~ x or discriminability ~ x and bias ~ x.")
      }
    } else {
      stop("Cannot interpret formula input. Please specify formulas as either
           discriminability = ~ x and bias = ~ x or discriminability ~ x and bias ~ x.")
    }
  } else form_bias_std <- form_bias


  return(list(form_disc_std, form_bias_std))

}

