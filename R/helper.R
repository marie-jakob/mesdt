.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  options("mesdt.backend" = "lme4")
  packageStartupMessage(
    sprintf("mesdt version %s\nSetting backend to lme4.", ver)
  )
}


#' Set backend for mesdt.
#'
#' @param backend_name `character` indicating the backend that should be used
#' to estimate the models (`lme4`, the default, which is set when loading the
#' package, or `glmmTMB`).
#'
#' @return NULL
#' @description
#' Allows the user to change the package that is used to estimate the
#' models. When loading the package, the backend is automatically set to `lme4`
#' (which is a dependency for this package). If `glmmTMB` is installed, the
#' backend can be changed to that package through this function.
#'
#' @export
#'
#' @examples
#' set_backend("lme4")
#' \dontrun{
#' set_backend("glmmTMB")
#' }
set_backend <- function(backend_name) {
  backend_name <- match.arg(backend_name, c("lme4", "glmmTMB"))
  if (! requireNamespace(backend_name, quietly = TRUE)) {
    msg <- paste("The", backend_name, "package must be installed to be used as
                 a backend.")
    stop(msg)
  } else {
    message(paste("Setting backend to ", backend_name, ".", sep = ""))
    options("mesdt.backend" = backend_name)
  }
}



#' Get the backend for mesdt.
#'
#'
#' @return the current backend used by mesdt.
#' @description
#'
#' Allows the user to check which package is currently set as the backend
#' to the estimate the models.
#'
#' @export
#'
#' @examples
#' get_backend()
get_backend <- function() {
    return(options("mesdt.backend"))
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


normalize_terms <- function(formula) {
  tl <- attr(terms(formula), "term.labels")

  # Split interaction terms (containing :) and sort components
  normalized <- sapply(tl, function(term) {
    if (grepl(":", term)) {
      parts <- unlist(strsplit(term, ":"))
      paste(sort(parts), collapse = ":")
    } else {
      term
    }
  })

  sort(normalized)  # Sort terms for comparison
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
  else if (x %in% c(3, "III", "iii", "3", "three")) type_std <- 3
  else type_std <- NULL
  return(type_std)
}


standardize_tests_input <- function(x) {
  if (x %in% c("LRT", "lrt", "LRTs")) test_std <- "lrt"
  else if (x %in% c("boot", "Boot", "bootstrap", "Bootstrap", "pb", "PB")) test_std <- "bootstrap"
  else test_std <- NULL
  return(test_std)
}



check_sensitivity <- function(fit_obj) {
  summ_mesdt <- summary(fit_obj)
  mu_mean <- summ_mesdt$d_coef[rownames(summ_mesdt$d_coef) == "(Intercept)", 1]
  if (length(mu_mean) > 0) {
    if (mu_mean < 0) warning("Mean population discriminability is < 0, indicating that the trial_type or response ('dv') variables might be coded reversely (see the documentation)")
  }
}

code_binary_var <- function(data, var, arg, role, noise_value, new_col, level = NULL) {
  if (!is.character(var) || length(var) != 1)
    stop("'", arg, "' must be a single character string naming a variable in data.")
  if (!var %in% names(data))
    stop("Given ", arg, " variable '", var, "' not in data.")
  x <- data[[var]]

  x_obs <- x[!is.na(x)]
  n_values <- length(unique(x_obs))

  if (n_values != 2)
    stop(arg, " must be a binary variable, but '", var, "' has ",
         n_values, " distinct value", if (n_values != 1) "s", ".")

  if (!is.null(level) && !is.factor(x))
    stop("'", arg, "_signal_level' can only be used if ", arg, " is a factor.")

  if (is.numeric(x)) {
    vals <- sort(unique(x_obs))
    if (!(all(vals == c(0, 1)) || all(vals == c(-1, 1))))
      stop("If ", arg, " is numeric, it must code ", role, " as 1 and the other value as 0 or -1 ",
           "(found: ", paste(vals, collapse = ", "), ").")
    signal_level <- 1
    is_signal <- x == 1
  } else if (is.factor(x)) {
    if (nlevels(x) != 2)
      stop(arg, " '", var, "' must have exactly two factor levels, but has ",
           nlevels(x), " (", paste(levels(x), collapse = ", "), "). ",
           "Use droplevels() and set the contrasts afterwards.")
    if (!is.null(level)) {
      if (!is.character(level) || length(level) != 1 || !level %in% levels(x))
        stop("'", arg, "_signal_level' must be one of the levels of '", var, "': ",
             paste(levels(x), collapse = ", "), ".")
      signal_level <- level
      src <- " (as specified)"
    } else {
      cm <- stats::contrasts(x)
      if (sum(cm[, 1] == 1) != 1 || !all(cm[, 1] %in% c(-1, 0, 1)))
        stop("Cannot tell which level of '", var, "' codes ", role, ": exactly one level ",
             "must be coded 1 in contrasts(", var, "), the other 0 or -1. Current coding: ",
             paste(levels(x), "=", round(cm[, 1], 2), collapse = ", "), ".")
      signal_level <- levels(x)[cm[, 1] == 1]
      src <- if (is.null(attr(x, "contrasts")))
        paste0(" (no contrasts set; global option: ", getOption("contrasts")[["unordered"]], ")") else ""
    }
    is_signal <- x == signal_level
    message(arg, " '", var, "': level '", signal_level, "' = ", role, src, ".")
  } else {
    stop(arg, " must be a binary numeric variable or a factor, but '", var,
         "' is of class '", class(x)[1], "'.")
  }
  data[[new_col]] <- ifelse(is_signal, 1, noise_value)
  return(list(data = data, var = new_col, signal_level = signal_level))
}
check_dv <- function(data, dv, dv_signal_level = NULL)
  code_binary_var(data, dv, arg = "dv", role = "signal responses",
                  noise_value = 0, new_col = "dv_num", level = dv_signal_level)

check_trial_type_var <- function(data, trial_type_var, trial_type_signal_level = NULL)
  code_binary_var(data, trial_type_var, arg = "trial_type", role = "signal trials",
                  noise_value = -1, new_col = "trial_type_num", level = trial_type_signal_level)


check_fit_formulas <- function(form_disc, form_bias) {
  # Check if formula has a left-hand side
  form_disc_std <- ""
  form_bias_std <- ""
  if (length(as.character(form_disc)) == 3) {
    stop("Cannot interpret formula input for discriminability. Please use
         one-sided formulas such as discriminability = ~ x.")
  } else if (length(as.character(form_bias)) == 3) {
    stop("Cannot interpret formula input for response bias. Please use
         one-sided formulas such as response_bias = ~ x.")
  }
  return()
}


#' @importFrom stats aggregate
#' @importFrom stats reshape
aggregate_data <- function(data, form_disc, form_bias, dv, trial_type) {
  grouping_vars <- c(unique(c(all.vars(form_disc), all.vars(form_bias))),
                     trial_type, dv)
  by_all <- data[, which(names(data) %in% grouping_vars )]
  dat_agg <- stats::aggregate(data[[dv]], by = by_all, FUN = length)
  id_vars <- c(unique(c(all.vars(form_disc), all.vars(form_bias))), trial_type)
  dat_agg_wide <- reshape(
    dat_agg,
    idvar = id_vars,          # column(s) identifying unique rows
    timevar = dv,    # column whose values become new columns
    direction = "wide"
  )
  new_names <- c("x.0", "x.1")

  dat_agg_wide[[new_names[1]]][is.na(dat_agg_wide[[new_names[1]]])] <- 0
  dat_agg_wide[[new_names[2]]][is.na(dat_agg_wide[[new_names[2]]])] <- 0


  dat_agg_wide$weight <- dat_agg_wide[[new_names[1]]] + dat_agg_wide[[new_names[2]]]
  dat_agg_wide$dv_agg <- dat_agg_wide[[new_names[2]]] / dat_agg_wide$weight
  return(dat_agg_wide)
}

#' @importFrom stats contrasts
#' @importFrom stats contr.sum
check_predictors <- function(data, form_disc, form_bias) {
  all_preds <- unique(all.vars(lme4::nobars(form_disc)), all.vars(lme4::nobars(form_bias)))
  not_sum_contrasts <- c()
  not_centered <- c()
  not_num_or_fac <- c()
  for (pred in all_preds) {
    if (is.numeric(data[[pred]])) {
      if (abs(mean(data[[pred]])) > 1e-5) {
        not_centered <- c(not_centered, pred)
      }
    } else if (inherits(data[[pred]], "factor")) {
      if (! all(sort(contrasts(data[[pred]])) == sort(contr.sum(nlevels(data[[pred]]))))) {
        not_sum_contrasts <- c(not_sum_contrasts, pred)
      }
    } else {
      not_num_or_fac <- c(not_num_or_fac, pred)
    }
  }
  if (length(not_sum_contrasts > 0)) {
    warning(
      paste0(
        "The following categorical predictors do not use sum contrasts: ",
        paste(not_sum_contrasts, collapse = ", "),
        "\nWe strongly recommend using sum contrasts for a straightforward interpretation ",
        "of lower-order effects in the presence of higher-order terms."
      )
    )
  }

  if (length(not_centered > 0)) {
    warning(
      paste0(
        "The following numerical predictors are not centered at zero: ",
        paste(not_centered, collapse = ", "),
        "\nWe strongly recommend centering predictors to simplify interpretation."
      )
    )
  }
  if (length(not_num_or_fac > 0)) {
    warning(
      paste0(
        "The following predictors are neither numeric nor factors: ",
        paste(not_num_or_fac, collapse = ", "),
        "\nPlease set the datatype of all predictors explicitly to prevent
        unexpected behavior of the fitting functions."
      )
    )
  }
  return()
}


handle_missing <- function(data, vars, na.action) {
  vars <- intersect(unique(vars), names(data))
  kept <- match.fun(na.action)(data[, vars, drop = FALSE])  # na.fail stops here
  omitted <- attr(kept, "na.action")
  if (!is.null(omitted)) data <- data[-omitted, , drop = FALSE]
  list(data = data, n_removed = length(omitted))
}

#' Remove trials with missing values in the variables used by the model
#' @keywords internal
remove_incomplete_trials <- function(data, vars) {
  vars <- intersect(unique(vars), names(data))
  na_counts <- colSums(is.na(data[, vars, drop = FALSE]))
  complete <- stats::complete.cases(data[, vars, drop = FALSE])
  if (!all(complete)) {
    message("Removed ", sum(!complete), " of ", nrow(data), " trials with missing values in: ",
            paste(names(na_counts)[na_counts > 0], collapse = ", "), ".")
    data <- data[complete, , drop = FALSE]
  }
  data
}

# TRUE for coefficient names belonging to the given SDT parameter ("mu" or "lambda")
is_sdt_par <- function(names, par) startsWith(names, paste0('mm[["', par, '"]]'))
