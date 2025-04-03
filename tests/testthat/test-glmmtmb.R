
options("mesdt.backend" = "lme4")

#------------------------------------------------------------------------------#
#### fit_mesdt() ####


test_that("fit_mesdt() estimates the correct model with glmmTMB", {

  options("mesdt.backend" = "glmmTMB")

  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)[[1]]), length(fixef(model_test_tmb)[[1]]))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)[[1]]), length(ranef(model_test_tmb)[[1]]))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test_tmb))))

  # dfs & LL
  expect_equal(df.residual(fit), df.residual(model_test_tmb))
  expect_equal(logLik(fit), logLik(model_test_tmb), tolerance = 1e-5)

  # fixed effects estimates
  expect_equal(as.numeric(unlist(fixef(fit))), as.numeric(unlist(fixef(model_test_tmb))), tolerance = 1e-5)

  # observed Fisher information
  expect_equal(unname(vcov(fit)[[1]][1:4, 1:4]), unname(vcov(model_test_tmb)[[1]][1:4, 1:4]), tolerance = 1e-5)

  # random effect variances and covariance
  expect_equal(unlist(VarCorr(fit)), unlist(VarCorr(model_test_tmb)), tolerance = 1e-4)

  # random effects estimates
  expect_equal(as.numeric(unlist(ranef(fit))), as.numeric(unlist(ranef(model_test_tmb))), tolerance = 1e-4)
}
)


# Results differ to some small extent (sometimes I need tolerances of 1e-1 for the tests to work)
# glmmTMB seems to get smaller likelihoods


test_that("glmmTMB and lme4 get similar results for fitted models and LRTs (crossed random effects)", {
  options("mesdt.backend" = "glmmTMB")
  form_mu <- ~ committee_ef + (1 | id) + (1 | file_name)
  form_lambda <- ~ committee_ef + (committee_ef | id)

  options("mesdt.backend" = "lme4")
  fit_lme <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts_lme <- compute_tests(fit_lme,
                                    data = dat_exp_2,
                                    type = 2,
                                    test_intercepts = T)

  # Type 2 - test_intercepts = F
  LRTs_2_lme <- compute_tests(fit_lme,
                         data = dat_exp_2,
                         type = 2,
                         test_intercepts = F)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts_lme <- compute_tests(fit_lme,
                                    data = dat_exp_2,
                                    type = 3,
                                    test_intercepts = T)

  # Type 3 - test_intercepts = F
  LRTs_3_lme <- compute_tests(fit_lme,
                              data = dat_exp_2,
                              type = 3,
                              test_intercepts = F)

  options("mesdt.backend" = "glmmTMB")
  fit_tmb <- fit_mesdt(form_mu, form_lambda,
                       dv = "assessment",
                       trial_type_var = "status_fac",
                       data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts_tmb <- compute_tests(fit_tmb,
                                    data = dat_exp_2,
                                    type = 2,
                                    test_intercepts = T)

  # Type 2 - test_intercepts = F
  LRTs_2_tmb <- compute_tests(fit_tmb,
                         data = dat_exp_2,
                         type = 2,
                         test_intercepts = F)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts_tmb <- compute_tests(fit_tmb,
                                         data = dat_exp_2,
                                         type = 3,
                                         test_intercepts = T)

  # Type 3 - test_intercepts = F
  LRTs_3_tmb <- compute_tests(fit_tmb,
                              data = dat_exp_2,
                              type = 3,
                              test_intercepts = F)


  expect_equal(fixef(fit_tmb$fit_obj)[[1]], fixef(fit_lme$fit_obj), tolerance = 1e-1)
  expect_equal(logLik(fit_tmb$fit), logLik(fit_lme$fit), tolerance = 1e-3)

  # LRTs
  expect_equal(chisquares_cross_2, as.numeric(LRTs_2_intercepts_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-3)
  expect_equal(chisquares_cross_2[c(2, 4)], as.numeric(LRTs_2_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)
  expect_equal(chisquares_cross_3, as.numeric(LRTs_3_intercepts_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)
  expect_equal(chisquares_cross_3[c(2, 4)], as.numeric(LRTs_3_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)

  expect_equal(as.numeric(LRTs_2_intercepts_lme$LRTs$LRT_results[, 4]),
               as.numeric(LRTs_2_intercepts_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-3)
  expect_equal(as.numeric(LRTs_2_lme$LRTs$LRT_results[, 4]),
               as.numeric(LRTs_2_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)
  expect_equal(as.numeric(LRTs_3_intercepts_lme$LRTs$LRT_results[, 4]),
               as.numeric(LRTs_3_intercepts_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)
  expect_equal(as.numeric(LRTs_3_lme$LRTs$LRT_results[, 4]),
               as.numeric(LRTs_3_tmb$LRTs$LRT_results[, 4]), tolerance = 1e-2)

})
