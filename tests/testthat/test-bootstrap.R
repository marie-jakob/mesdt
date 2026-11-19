
set_backend("glmmTMB")

test_that("compute_tests() compares the correct models for bootstrap tests of fixed effects", {
    fit <- fit_mesdt(~ committee * emp_gender + (1 | id),
                     ~ committee * emp_gender_ef + (1 | id),
                     dv = "assessment",
                     trial_type = "status_fac",
                     data = dat_exp_2)
        suppressWarnings(boot <- compute_tests(fit,
                              tests = "bootstrap",
                              type = 2,
                              nsim = 1,
                              test_intercepts = TRUE))
        expect_equal(boot$pb_test_results[, 1], boot$LRT_results[, 4])
  })


test_that("compute_tests() compares the correct models for bootstrap tests of factors with multiple levels", {

    form_mu <- ~ contingencies + (1 | id)
    form_lambda <- ~ contingencies + (1 | id)
    # Same for the uncorrelated model
    fit <- fit_mesdt(form_mu, form_lambda,
                     dv = "assessment",
                     trial_type = "status_fac",
                     data = dat_exp_2)
        suppressWarnings(boot <- compute_tests(fit,
                              test_intercepts = FALSE,
                              type = 3,
                              tests = "bootstrap",
                              nsim = 1))
        expect_equal(boot$pb_test_results[, 1], boot$LRT_results[, 4])

})

test_that("parametric bootstrap works for aggregated data with glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  p_boot <- function(backend) {
    set_backend(backend)
    fit <- fit_mesdt(~ 1 + (1 | id), ~ committee + (1 | id), dv = "assessment",
                     data = debi3subset, trial_type = "status", aggregate = TRUE)
    compute_tests(fit, tests = "bootstrap", nsim = 30, seed = 1,
                  tests_discriminability = NULL)$pb_test_results[, "p.value"]
  }
  p_lme4 <- p_boot("lme4")
  p_tmb <- p_boot("glmmTMB")
  set_backend("lme4")
  expect_equal(as.numeric(p_tmb), as.numeric(p_lme4), tolerance = 0.1)
})

