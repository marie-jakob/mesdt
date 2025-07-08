
options("mesdt.backend" = "lme4")

#------------------------------------------------------------------------------#
#### compute_tests() ####

test_that("compute_tests() computes the correct Chisq value for correlated random effects.", {
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type = "trial_type_fac")
  lrts_test <- compute_tests(fit, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRT_results[, 4])), model_test_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRT_results[, 5])), model_test_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

test_that("compute_tests() computes the correct Chisq value for uncorrelated random effects.", {
  # Same for the uncorrelated model
  fit <- fit_mesdt(~ x1 + (x1 || ID), ~ x1 + (x1 || ID), dv = "y",
                   trial_type = "trial_type_fac", data = internal_sdt_data)
  lrts_test <- compute_tests(fit, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRT_results[, 4])), model_test_uncor_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRT_results[, 5])), model_test_uncor_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

#------------------------------------------------------------------------------#
#### Type II SS ####


test_that("compute_tests() works only with intercepts", {

  form_mu <- ~ 1 + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  fit <- fit_mesdt(form_lambda, form_mu,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = T)
  LRTs_3 <- compute_tests(fit,
                         type = 3,
                         test_intercepts = T)
  # Should be equal to each other
  expect_equal(unlist(LRTs_3$LRT_results[, 4]), unlist(LRTs_2$LRT_results[, 4]), tolerance = 1e-4)
  # And equal to models fitted by hand
  expect_equal(as.numeric(LRTs_2$LRT_results[, 4]), chi_squares_intercepts, tolerance = 1e-4)
})



test_that("compute_tests() Type II works with one predictor on mu", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                       type = 2,
                       test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRT_results[, 4]), chi_squares_one_pred_mu_2, tolerance = 1e-4)

})


##############################
test_that("compute_tests() Type II works with one predictor on mu and lambda", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ committee_ef + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)


  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRT_results[, 4]), chisquares_one_factor_2[c(2, 4)], tolerance = 1e-4)
})


test_that("compute_tests() works with a standard two-factorial design", {
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias =~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)


  # Type III, test_intercepts = T
  fit_LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = T)
  expect_equal(chisquares_two_factors_3, as.numeric(fit_LRTs$LRT_results[, 4]), tolerance = 1e-3)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[7]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$reduced_fits[[8]]), stats::df.residual(fit$fit_obj) + 1)

})



#------------------------------------------------------------------------------#
#### Factors with > 2 levels ####

test_that("compute_tests() works for factors with > 2 levels", {
  # Type II, test_intercepts = T
  form_mu <- ~ contingencies + (1 | id)
  form_lambda <- ~ contingencies + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)
  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRT_results[, 4]), chisquares_contingencies_3[c(2, 4)], tolerance = 1e-4)
})

#------------------------------------------------------------------------------#
#### Crossed random effects ####

test_that("compute_tests() Type II works with one predictor on mu and lambda and crossed random effects", {
  form_mu <- ~ committee_ef + (1 | id) + (1 | file_name)
  form_lambda <- ~ committee_ef + (committee_ef | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                                    type = 2,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRT_results[, 4]), chisquares_cross_2, tolerance = 1e-3)
})


#------------------------------------------------------------------------------#
#### Only test selected random effects ####

test_that("compute_tests() works for only selected effects", {
  options("mesdt.backend" = "lme4")
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                       type = 2,
                       test_intercepts = T,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)

  expect_equal(chisquares_two_factors_2[c(1, 2, 5, 7)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 2,
                       test_intercepts = F,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)
  expect_equal(chisquares_two_factors_2[c(2, 7)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = T,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)

  expect_equal(chisquares_two_factors_3[c(1, 2, 5, 7)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = F,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)
  expect_equal(chisquares_two_factors_3[c(2, 7)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

})


#------------------------------------------------------------------------------#
#### compute_tests() works only for mu and lambda ####

test_that("compute_tests() works for only tests on mu", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_discriminability = ~ committee,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_2[c(5, 6)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # Type III
  LRTs <- compute_tests(fit,
                        type = 3,
                        test_intercepts = T,
                        tests_discriminability = ~ committee,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_3[c(5, 6)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)
})


test_that("compute_tests() works for only tests on lambda", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_response_bias = ~ committee,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[c(1, 2)], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # all parameters on mu, none on lambda
  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[1:4], as.numeric(LRTs$LRT_results[, 4]), tolerance = 1e-4)
})



#------------------------------------------------------------------------------#
#### Control Arguments on compute_tests() ####

test_that("compute_tests() computes the correct Chisq value for correlated random effects.", {
  skip_if_not_installed("glmmTMB")
  library(glmmTMB)
  options("mesdt.backend" = "lme4")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type = "trial_type_fac",
                   control = lme4::glmerControl(optCtrl = list(maxfun = 1234)))
  lrts_test <- compute_tests(fit, test_intercepts = F,
                             control = lme4::glmerControl(optCtrl = list(maxfun = 1234)))
  expect_equal(lrts_test$reduced_fits$x1_lambda@optinfo$control$maxfun, 1234)

  options("mesdt.backend" = "glmmTMB")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type = "trial_type_fac",
                   control = glmmTMB::glmmTMBControl(list(iter.max=12334, eval.max=12334)))
  lrts_test <- compute_tests(fit, test_intercepts = F,
                             control = glmmTMB::glmmTMBControl(list(iter.max=12334, eval.max=12334)))
  expect_true("control" %in% as.character(lrts_test$reduced_fits$x1_lambda$call))

})

options("mesdt.backend" = "lme4")


#------------------------------------------------------------------------------#
#### Check that error message is given when the formula contain terms not in the model ####

test_that("compute_tests() throws an error when the formula contain terms not in the model", {
  options("mesdt.backend" = "lme4")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type = "trial_type_fac")
  expect_error(lrts_test <- compute_tests(fit, tests_discriminability = ~ x2,
                                          test_intercepts = T))
  expect_error(lrts_test <- compute_tests(fit, tests_response_bias = ~ x2,
                                          test_intercepts = T))
})


