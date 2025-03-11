.onLoad <- function(libname, pkgname) {
  options(mlsdt.backend = "lme4")
}
fit_glmm <- function(glmer_formula,
                     data,
                     mm,
                     control = NULL) {
  # mm <- construct_modelmatrices(formula_mu, formula_lambda, dv, data, trial_type_var)
  # get global options
  if (! (options("mlsdt.backend") %in% c("lme4", "glmmTMB"))) {
    message("Only lme4 and glmmTMB backends are supported at the moment. Defaulting to lme4.")
    (options("mlsdt.backend" = "lme4"))
  }

  if (options("mlsdt.backend") == "glmmTMB" & !requireNamespace("glmmTMB", quietly = TRUE)) {
    message("Package \"glmmTMB\" must be installed to use it as backend. Setting backend to lme4.")
    options("mlsdt.backend" = "lme4")
  }

  # only applies to lme4, defaults to 0 at the moment
  # -> might be changed later to a function argument
  # nAGQ <<- ifelse(is.null(options("nAGQ")$nAGQ), 0, options("nAGQ"))
  #print(nAGQ)


  # Fit a GLM if there are no random effects
  if (is.null(lme4::findbars(glmer_formula))) {
    message("Formula does not contain any random terms. Fitting a single-level model.")
    fit_obj <- stats::glm(glmer_formula,
                          data = data,
                          family = binomial(link = "probit"))
  } else if ((options("mlsdt.backend") == "lme4")) {
    if (is.null(control) |
        ! "glmerControl" %in% attr(control, "class")) {
      if (! is.null(control)) message("Invalid control argument for lme4. Using the default settings.")
      fit_obj <- lme4::glmer(glmer_formula,
                             data = data,
                             family = binomial(link = "probit"),
                             # this is only for testing speed -> changed for actual use
                             nAGQ = 0)
    } else {
      fit_obj <- lme4::glmer(glmer_formula,
                             data = data,
                             family = binomial(link = "probit"),
                             # this is only for testing speed -> changed for actual use
                             nAGQ = 0,
                             control = control)
    }

  } else if ((options("mlsdt.backend") == "glmmTMB")) {
    #print("fitting with glmmTMB")
    if (is.null(control) |
        "glmerControl" %in% attr(control, "class")) {
      if (! is.null(control)) message("Invalid control argument for glmmTMB. Using the default settings.")
      fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                  data = data,
                                family = binomial(link = "probit"))
    } else {
      fit_obj <- glmmTMB::glmmTMB(glmer_formula,
                                  data = data,
                                  family = binomial(link = "probit"),
                                  control = control)
    }
  }
  return(fit_obj)
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
  print(matches)
  print(matches > 0)
  return(matches[matches > 0])
}

