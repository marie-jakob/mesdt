
options("mesdt.backend" = "lme4")

#------------------------------------------------------------------------------#
#### compute_tests() ####

test_that("compute_tests() computes the correct Chisq value for correlated random effects.", {
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac")
  lrts_test <- compute_tests(fit, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 4])), model_test_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 5])), model_test_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

test_that("compute_tests() computes the correct Chisq value for uncorrelated random effects.", {
  # Same for the uncorrelated model
  fit <- fit_mesdt(~ x1 + (x1 || ID), ~ x1 + (x1 || ID), dv = "y",
                   trial_type_var = "trial_type_fac", data = internal_sdt_data)
  lrts_test <- compute_tests(fit, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 4])), model_test_uncor_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 5])), model_test_uncor_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

test_that("compute_tests() throws a message when there is nothing to test in the model.", {
  # Case 1: no predictors & test_intercepts = F
  fit <- fit_mesdt(~ 1 + (x1 || ID), ~ 1 + (x1 || ID), dv = "y",
                   trial_type_var = "trial_type_fac", data = internal_sdt_data)
  expect_message(compute_tests(fit, test_intercepts = F))

  # TODO: this does not throw the expected message
})


#------------------------------------------------------------------------------#
#### Type II SS ####


test_that("compute_tests() works only with intercepts", {

  form_mu <- ~ 1 + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  fit <- fit_mesdt(form_lambda, form_mu,
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = T)
  LRTs_3 <- compute_tests(fit,
                         type = 3,
                         test_intercepts = T)
  # Should be equal to each other
  expect_equal(unlist(LRTs_3$LRTs$LRT_results[, 4]), unlist(LRTs_2$LRTs$LRT_results[, 4]), tolerance = 1e-4)
  # And equal to models fitted by hand
  expect_equal(as.numeric(LRTs_2$LRTs$LRT_results[, 4]), chi_squares_intercepts, tolerance = 1e-4)
})



test_that("compute_tests() Type II works with one predictor on mu", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                       type = 2,
                       test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs$LRT_results[, 4]), chi_squares_one_pred_mu_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_tests(fit,
                         type = 3,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs$LRT_results[, 4]), chi_squares_one_pred_mu_2[3], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chi_squares_one_pred_mu_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chi_squares_one_pred_mu_3[3], tolerance = 1e-4)

})


##############################
test_that("compute_tests() Type II works with one predictor on mu and lambda", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ committee_ef + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                                    type = 2,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs$LRT_results[, 4]), chisquares_one_factor_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs$LRT_results[, 4]), chisquares_one_factor_2[c(2, 4)], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chisquares_one_factor_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chisquares_one_factor_3[c(2, 4)], tolerance = 1e-4)

})


test_that("compute_tests() works with a standard two-factorial design", {
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias =~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  fit_LRTs <- compute_tests(fit,
                       type = 2,
                       test_intercepts = T)
  expect_equal(chisquares_two_factors_2, as.numeric(fit_LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)
  # check dfs
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 4)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 4)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[7]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[8]]), stats::df.residual(fit$fit_obj) + 1)

  # Type II, test_intercepts = F
  fit_LRTs <- compute_tests(fit,
                       type = 2)
  expect_equal(chisquares_two_factors_2[-c(1, 5)], as.numeric(fit_LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)


  # Type III, test_intercepts = T
  fit_LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = T)
  expect_equal(chisquares_two_factors_3, as.numeric(fit_LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-3)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[7]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[8]]), stats::df.residual(fit$fit_obj) + 1)

  # Type III, test_intercepts = T
  fit_LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = F)
  expect_equal(chisquares_two_factors_3[-c(1, 5)], as.numeric(fit_LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-3)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(fit_LRTs$LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)
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
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                                    type = 2,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs$LRT_results[, 4]), chisquares_contingencies_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs$LRT_results[, 4]), chisquares_contingencies_2[c(2, 4)], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chisquares_contingencies_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chisquares_contingencies_3[c(2, 4)], tolerance = 1e-4)

})

#------------------------------------------------------------------------------#
#### Crossed random effects ####

test_that("compute_tests() Type II works with one predictor on mu and lambda and crossed random effects", {
  form_mu <- ~ committee_ef + (1 | id) + (1 | file_name)
  form_lambda <- ~ committee_ef + (committee_ef | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_tests(fit,
                                    type = 2,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs$LRT_results[, 4]), chisquares_cross_2, tolerance = 1e-3)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_tests(fit,
                         type = 2,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs$LRT_results[, 4]), chisquares_cross_2[c(2, 4)], tolerance = 1e-3)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_tests(fit,
                                    type = 3,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs$LRT_results[, 4]), chisquares_cross_3, tolerance = 1e-3)

  # Type 3 - test_intercepts = F
  LRTs_3 <- compute_tests(fit,
                          type = 3,
                          test_intercepts = F)
  expect_equal(as.numeric(LRTs_3$LRTs$LRT_results[, 4]), chisquares_cross_3[c(2, 4)], tolerance = 1e-3)

})


#------------------------------------------------------------------------------#
#### Random Effects ####

test_that("compute_tests() works for testing random effects", {
  options("mesdt.backend" = "glmmTMB")
  # with correlations
  fit <- fit_mesdt(~ committee + (1 | id), ~ committee + (1 | id), dv = "assessment", data = dat_exp_2,
                   correlate_sdt_params = T,
                   trial_type_var = "status_fac")

  lrts_test <- compute_tests(fit, test_intercepts = T, test_ran_ef = T)

  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 4]), chi_squares_rdm_intercepts)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 5]), pchisqmix(chi_squares_rdm_intercepts, 2, 0.5, lower.tail = F))

  # without correlations
  fit <- fit_mesdt(~ committee + (1 | id),
                   ~ committee + (1 | id),
                   dv = "assessment",
                   data = dat_exp_2,
                   correlate_sdt_params = F,
                   trial_type_var = "status_fac")

  lrts_test <- compute_tests(fit, test_intercepts = T, test_ran_ef = T)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 4]), chi_squares_rdm_intercepts_uncor)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 5]), pchisqmix(chi_squares_rdm_intercepts_uncor, 1, 0.5, lower.tail = F))
})


# TODO: does the random stuff work for factors with multiple levels?

test_that("compute_tests() works for testing crossed random effects", {
  options("mesdt.backend" = "glmmTMB")
  # with correlations
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                   trial_type_var = "status_fac")
  lrts_test <- compute_tests(fit, test_intercepts = T, test_ran_ef = T)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 4]), chi_squares_rdm_cross, tolerance = 1e-5)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[1:3, 5]), pchisqmix(chi_squares_rdm_cross[1:3], 3, 0.5, lower.tail = F))
  expect_equal(as.numeric(lrts_test$LRTs$LRT_results[4, 5]), pchisqmix(chi_squares_rdm_cross[4], 1, 0.5, lower.tail = F))

  # without correlations
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee || id), dv = "assessment", data = dat_exp_2,
                   trial_type_var = "status_fac")
  lrts_test <- compute_tests(fit, test_intercepts = T, test_ran_ef = T)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[, 4]), chi_squares_rdm_cross_uncor, tolerance = 1e-5)
  expect_equal(unlist(lrts_test$LRTs$LRT_results[1:3, 5]), pchisqmix(chi_squares_rdm_cross_uncor[1:3], 1, 0.5, lower.tail = F))
  expect_equal(as.numeric(lrts_test$LRTs$LRT_results[4, 5]), pchisqmix(chi_squares_rdm_cross_uncor[4], 1, 0.5, lower.tail = F))

})


test_that("compute_tests() works for testing crossed random effects without the intercept", {
  options("mesdt.backend" = "glmmTMB")
  # with correlations
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,trial_type_var = "status_fac")
  lrts_test <- compute_tests(fit, test_intercepts = F, test_ran_ef = T)
  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 4])), chi_squares_rdm_cross[2], tolerance = 1e-5)

  # without correlations
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee || id), dv = "assessment", data = dat_exp_2,
                   trial_type_var = "status_fac")
  lrts_test <- compute_tests(fit, test_intercepts = F, test_ran_ef = T)
  expect_equal(unname(unlist(lrts_test$LRTs$LRT_results[, 4])), chi_squares_rdm_cross_uncor[2], tolerance = 1e-5)

})


#------------------------------------------------------------------------------#
#### Only test selected random effects ####

test_that("compute_tests() works for only selected effects", {
  options("mesdt.backend" = "lme4")
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                       type = 2,
                       test_intercepts = T,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)

  expect_equal(chisquares_two_factors_2[c(1, 2, 5, 7)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 2,
                       test_intercepts = F,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)
  expect_equal(chisquares_two_factors_2[c(2, 7)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = T,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)

  expect_equal(chisquares_two_factors_3[c(1, 2, 5, 7)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  LRTs <- compute_tests(fit,
                       type = 3,
                       test_intercepts = F,
                       tests_response_bias = ~ committee,
                       tests_discriminability = ~ emp_gender)
  expect_equal(chisquares_two_factors_3[c(2, 7)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

})

#------------------------------------------------------------------------------#
#### compute_tests() works only for mu and lambda ####

test_that("compute_tests() works for only tests on mu", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_discriminability = ~ committee,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_2[c(5, 6)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # Type III
  LRTs <- compute_tests(fit,
                        type = 3,
                        test_intercepts = T,
                        tests_discriminability = ~ committee,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_3[c(5, 6)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # No intercept
  LRTs <- compute_tests(fit,
                        test_intercepts = F,
                        tests_response_bias = ~ committee:emp_gender,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_3[c(4)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # all parameters on mu, none on lambda
  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_2[5:8], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)
})


test_that("compute_tests() works for only tests on lambda", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_response_bias = ~ committee,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[c(1, 2)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # Type III
  LRTs <- compute_tests(fit,
                        type = 3,
                        test_intercepts = T,
                        tests_response_bias = ~ committee,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_3[c(1, 2)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # all parameters on mu, none on lambda
  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[1:4], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)
})


#------------------------------------------------------------------------------#
#### compute_tests() works multiple parameters (but not all) ####

test_that("compute_tests() works on multiple parameters", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  expect_warning(LRTs <- compute_tests(fit,
                               type = 2,
                               test_intercepts = T,
                               tests_discriminability = ~ committee + emp_gender + committee:emp_gender,
                               tests_response_bias = ~ committee + emp_gender + committee:emp_gender), regexp = NA)

  expect_warning(LRTs <- compute_tests(fit,
                                       type = 2,
                                       test_intercepts = T,
                                       tests_response_bias = NULL,
                                       tests_discriminability = ~ committee + emp_gender), regexp = NA)
  # Type III
  LRTs <- compute_tests(fit,
                        type = 3,
                        test_intercepts = T,
                        tests_discriminability = ~ committee,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_3[c(5, 6)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # No intercept
  LRTs <- compute_tests(fit,
                        test_intercepts = F,
                        tests_response_bias = ~ committee:emp_gender,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_3[c(4)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # all parameters on mu, none on lambda
  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_response_bias = NULL)

  expect_equal(chisquares_two_factors_2[5:8], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)
})

#------------------------------------------------------------------------------#
#### Only tests on lambda or mu ####

test_that("compute_tests() works for only tests on lambda", {
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)

  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_response_bias = ~ committee,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[c(1, 2)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # Type III
  LRTs <- compute_tests(fit,
                        type = 3,
                        test_intercepts = T,
                        tests_response_bias = ~ committee,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_3[c(1, 2)], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)

  # all parameters on mu, none on lambda
  LRTs <- compute_tests(fit,
                        type = 2,
                        test_intercepts = T,
                        tests_discriminability = NULL)

  expect_equal(chisquares_two_factors_2[1:4], as.numeric(LRTs$LRTs$LRT_results[, 4]), tolerance = 1e-4)
})



#------------------------------------------------------------------------------#
#### Control Arguments on compute_tests() ####

test_that("compute_tests() computes the correct Chisq value for correlated random effects.", {
  options("mesdt.backend" = "lme4")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac",
                   control = lme4::glmerControl(optCtrl = list(maxfun = 1234)))
  lrts_test <- compute_tests(fit, test_intercepts = T,
                             control = lme4::glmerControl(optCtrl = list(maxfun = 1234)))
  expect_equal(lrts_test$LRTs$reduced_fits$Intercept@optinfo$control$maxfun, 1234)

  options("mesdt.backend" = "glmmTMB")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac",
                   control = glmmTMB::glmmTMBControl(list(iter.max=12334, eval.max=12334)))
  lrts_test <- compute_tests(fit, test_intercepts = T,
                             control = glmmTMB::glmmTMBControl(list(iter.max=12334, eval.max=12334)))
  expect_true("control" %in% as.character(lrts_test$LRTs$reduced_fits$Intercept$call))

})

options("mesdt.backend" = "lme4")


#------------------------------------------------------------------------------#
#### Check that error message is given when the formula contain terms not in the model ####

test_that("compute_tests() throws an error when the formula contain terms not in the model", {
  options("mesdt.backend" = "lme4")
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac")
  expect_error(lrts_test <- compute_tests(fit, tests_discriminability = ~ x2,
                                          test_intercepts = T))
  expect_error(lrts_test <- compute_tests(fit, tests_response_bias = ~ x2,
                                          test_intercepts = T))
})


