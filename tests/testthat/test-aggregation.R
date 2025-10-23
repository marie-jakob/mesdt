test_that("aggregation gives the same results as usual", {
  for (backend in c("lme4", "glmmTMB")) {
    skip_if_not_installed(backend)
    set_backend(backend)
    fit <- fit_mesdt(~ x1 + (1 | ID), ~ x1 + (1 | ID), dv = "y", data = internal_sdt_data,
                     trial_type = "trial_type_fac",
                     aggregate = FALSE)
    s <- summary(fit)

    fit_agg <- fit_mesdt(~ x1 + (1 | ID), ~ x1 + (1 | ID), dv = "y", data = internal_sdt_data,
                         trial_type = "trial_type_fac",
                         aggregate = TRUE)
    s_agg <- summary(fit_agg)

    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-5)
    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-5)

    lrts <- compute_tests(fit, test_intercepts = T)
    lrts_agg <- compute_tests(fit_agg, test_intercepts = T)

    expect_equal(lrts$LRT_results[, c(3, 4)], lrts_agg$LRT_results[, c(3, 4)])
  }

})


#------------------------------------------------------------------------------#
#### Different formulas for response bias and sensitivity ####



test_that("aggregation gives the same results as usual", {
  for (backend in c("lme4", "glmmTMB")) {
    skip_if_not_installed(backend)
    set_backend(backend)
    fit <- fit_mesdt(~ committee + (1 | id), ~ emp_gender + (1 | id),
                     dv = "assessment", data = debi3subset,
                     trial_type = "status",
                     aggregate = FALSE)
    s <- summary(fit)

    fit_agg <- fit_mesdt(~ committee + (1 | id), ~ emp_gender + (1 | id),
                         dv = "assessment", data = debi3subset,
                         trial_type = "status",
                         aggregate = TRUE)
    s_agg <- summary(fit_agg)

    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-4)
    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-4)

    lrts <- compute_tests(fit, test_intercepts = T)
    lrts_agg <- compute_tests(fit_agg, test_intercepts = T)

    expect_equal(lrts$LRT_results[, c(3, 4)], lrts_agg$LRT_results[, c(3, 4)], tolerance = 1e-2)
  }
})



test_that("effects only in random effects", {
  for (backend in c("lme4", "glmmTMB")) {
    skip_if_not_installed(backend)
    set_backend(backend)
    fit <- fit_mesdt(~ 1 + (1 | id), ~ 1 + (committee | id),
                     dv = "assessment", data = debi3subset,
                     trial_type = "status",
                     aggregate = FALSE)
    s <- summary(fit)

    fit_agg <- fit_mesdt(~ 1 + (1 | id), ~ 1 + (committee | id),
                         dv = "assessment", data = debi3subset,
                         trial_type = "status",
                         aggregate = TRUE)
    s_agg <- summary(fit_agg)

    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-3)
    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-3)

    lrts <- compute_tests(fit, test_intercepts = T)
    lrts_agg <- compute_tests(fit_agg, test_intercepts = T)

    expect_equal(lrts$LRT_results[, c(3, 4)], lrts_agg$LRT_results[, c(3, 4)], tolerance = 1e-3)
  }
})



#------------------------------------------------------------------------------#
#### Crossed random effects ####


test_that("crossed random effects", {
  for (backend in c("lme4", "glmmTMB")) {
    skip_if_not_installed(backend)
    set_backend(backend)
    fit <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id),
                     dv = "assessment", data = debi3subset,
                     trial_type = "status",
                     aggregate = FALSE)
    s <- summary(fit)

    fit_agg <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id),
                         dv = "assessment", data = debi3subset,
                         trial_type = "status",
                         aggregate = TRUE)
    s_agg <- summary(fit_agg)

    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-3)
    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-3)

    lrts <- compute_tests(fit, test_intercepts = T)
    lrts_agg <- compute_tests(fit_agg, test_intercepts = T)

    expect_equal(lrts$LRT_results[, c(3, 4)], lrts_agg$LRT_results[, c(3, 4)], tolerance = 1e-3)
  }
})


#------------------------------------------------------------------------------#
#### Gumbel ####


test_that("aggregation works for gumbel-min", {
  for (backend in c("lme4", "glmmTMB")) {
    skip_if_not_installed(backend)
    set_backend(backend)
    fit <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id),
                     dv = "assessment", data = debi3subset,
                     trial_type = "status",
                     aggregate = FALSE,
                     distribution = "gumbel-min")
    s <- summary(fit)

    fit_agg <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id),
                         dv = "assessment", data = debi3subset,
                         trial_type = "status",
                         aggregate = TRUE,
                         distribution = "gumbel-min")
    s_agg <- summary(fit_agg)

    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-3)
    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-3)

    lrts <- compute_tests(fit, test_intercepts = T)
    lrts_agg <- compute_tests(fit_agg, test_intercepts = T)

    expect_equal(lrts$LRT_results[, c(3, 4)], lrts_agg$LRT_results[, c(3, 4)], tolerance = 1e-3)
  }
})


#------------------------------------------------------------------------------#
#### Test bootstrapping ####


#test_that("bootstrapping with aggregated data works", {
#  for (backend in c("lme4", "glmmTMB")) {
#    skip_if_not_installed(backend)
#    set_backend(backend)
#    fit <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id),
#                     dv = "assessment", data = debi3subset,
#                     trial_type = "status",
#                     aggregate = FALSE)
#    s <- summary(fit)
#
#    fit_agg <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name), ~ 1 + (committee | id#),
#                         dv = "assessment", data = debi3subset,
#                         trial_type = "status",
#                         aggregate = TRUE)
#    s_agg <- summary(fit_agg)
#
#    expect_equal(s$d_coef, s_agg$d_coef, tolerance = 1e-3)
#    expect_equal(s$c_coef, s_agg$c_coef, tolerance = 1e-3)
#    library(parallel)
#    cl <- makeCluster(8)
#
#    boot <- compute_tests(fit, tests = "bootstrap", test_intercepts = T,
#                          nsim = 500, seed = 1, cl = cl)
#    boot_agg <- compute_tests(fit_agg, tests = "bootstrap", test_intercepts = T,
#                              nsim = 500, seed = 1, cl = cl)
#
#    expect_equal(boot$LRT_results[, c(3, 4)], boot_agg$LRT_results[, c(3, 4)], #tolerance = 1e-3)
#  }
#})
