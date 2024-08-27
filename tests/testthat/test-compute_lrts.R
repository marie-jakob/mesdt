
#------------------------------------------------------------------------------#
#### compute_LRTs() ####

test_that("compute_LRTs() computes the correct Chisq value for correlated random effects.", {
  fit <- fit_mlsdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data)$fit_obj
  mm <- construct_modelmatrices(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), data = internal_sdt_data)
  lrts_test <- compute_LRTs(fit, ~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                            mm  = mm, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRTs[, 4])), model_test_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRTs[, 5])), model_test_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

test_that("compute_LRTs() computes the correct Chisq value for uncorrelated random effects.", {
  # Same for the uncorrelated model
  fit <- fit_mlsdt(~ x1 + (x1 || ID), ~ x1 + (x1 || ID), dv = "y", data = internal_sdt_data)$fit_obj
  mm <- construct_modelmatrices(~ x1 + (x1 || ID), ~ x1 + (x1 || ID), data = internal_sdt_data)
  lrts_test <- compute_LRTs(fit, ~ x1 + (x1 || ID), ~ x1 + (x1 || ID), dv = "y", data = internal_sdt_data,
                            mm  = mm, test_intercepts = T)

  # Chisq values
  # low tolerance because the package function fits with nAGQ = 0 (afex with nAGQ = 1)
  expect_equal(unname(unlist(lrts_test$LRTs[, 4])), model_test_uncor_afex$anova_table$Chisq, tolerance = 1e-2)

  expect_equal(unname(unlist(lrts_test$LRTs[, 5])), model_test_uncor_afex$anova_table$`Pr(>Chisq)`, tolerance = 1e-2)
})

test_that("compute_LRTs() throws a message when there is nothing to test in the model.", {
  # Case 1: no predictors & test_intercepts = F
  fit <- fit_mlsdt(~ 1 + (x1 || ID), ~ 1 + (x1 || ID), dv = "y", data = internal_sdt_data)$fit_obj
  mm <- construct_modelmatrices(~ 1 + (x1 || ID), ~ 1 + (x1 || ID), data = internal_sdt_data)
  expect_message(compute_LRTs(fit, ~ 1 + (x1 || ID), ~ 1 + (x1 || ID), dv = "y", data = internal_sdt_data,
                              mm  = mm, test_intercepts = F))
  #expect_equal(compute_LRTs(fit, ~ 1 + (x1 || ID), ~ 0 + (x1 || ID), dv = "y", data = internal_sdt_data,
  #                          mm  = mm, test_intercepts = T),
  #             NULL)

})


#------------------------------------------------------------------------------#
#### Type II SS ####


test_that("compute_LRTs() works only with intercepts", {

  form_mu <- ~ 1 + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  fit <- fit_mlsdt(form_lambda, form_mu,
                   dv = "assessment",
                   trial_type_var = "status_ef",
                   data = dat_exp_2)

  mm <- construct_modelmatrices(form_mu, form_lambda,
                                dv = "assessment",
                                trial_type_var = "status",
                                data = dat_exp_2)

  LRTs_2 <- compute_LRTs(fit$fit_obj,
                         form_mu, form_lambda,
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 2,
                         mm,
                         test_intercepts = T)
  LRTs_3 <- compute_LRTs(fit$fit_obj,
                         form_mu, form_lambda,
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 3,
                         mm,
                         test_intercepts = T)
  # Should be equal to each other
  expect_equal(unlist(LRTs_3$LRTs[, 4]), unlist(LRTs_2$LRTs[, 4]), tolerance = 1e-4)
  # And equal to models fitted by hand
  expect_equal(as.numeric(LRTs_2$LRTs[, 4]), chi_squares_intercepts, tolerance = 1e-4)
})



test_that("compute_LRTs() Type II works with one predictor on mu", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ 1 + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mlsdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_ef",
                   data = dat_exp_2)

  mm <- construct_modelmatrices(form_mu, form_lambda,
                                dv = "assessment",
                                trial_type_var = "status",
                                data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                       dv = "assessment",
                       data = dat_exp_2,
                       type = 2,
                       mm,
                       test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs[, 4]), chi_squares_one_pred_mu_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_LRTs(fit$fit_obj,
                         form_mu, form_lambda,
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 2,
                         mm,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs[, 4]), chi_squares_one_pred_mu_2[3], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chi_squares_one_pred_mu_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chi_squares_one_pred_mu_3[3], tolerance = 1e-4)

})



test_that("compute_LRTs() Type II works with one predictor on mu and lambda", {
  form_mu <- ~ committee_ef + (1 | id)
  form_lambda <- ~ committee_ef + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mlsdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_ef",
                   data = dat_exp_2)

  mm <- construct_modelmatrices(form_mu, form_lambda,
                                dv = "assessment",
                                trial_type_var = "status",
                                data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 2,
                                    mm,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs[, 4]), chisquares_one_factor_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_LRTs(fit$fit_obj,
                         form_mu, form_lambda,
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 2,
                         mm,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs[, 4]), chisquares_one_factor_2[c(2, 4)], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chisquares_one_factor_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chisquares_one_factor_3[c(2, 4)], tolerance = 1e-4)

})


test_that("compute_LRTs() works with a standard two-factorial design", {
  # Type II, test_intercepts = T
  fit <- fit_mlsdt(formula_lambda = ~ committee_ef * emp_gender_ef + (1 | id),
                   formula_mu = ~ committee_ef * emp_gender_ef + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_ef",
                   data = dat_exp_2)

  mm <- construct_modelmatrices(~ committee * emp_gender +  (1 | id),
                                ~ committee * emp_gender + (1 | id),
                                dv = "assessment",
                                trial_type_var = "status",
                                data = dat_exp_2)

  LRTs <- compute_LRTs(fit$fit_obj,
                       ~ committee * emp_gender + (1 | id),
                       ~ committee * emp_gender + (1 | id),
                       dv = "assessment",
                       data = dat_exp_2,
                       type = 2,
                       mm,
                       test_intercepts = T)
  expect_equal(chisquares_two_factors_2, as.numeric(LRTs$LRTs[, 4]), tolerance = 1e-4)
  # check dfs
  expect_equal(stats::df.residual(LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 4)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 4)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[7]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[8]]), stats::df.residual(fit$fit_obj) + 1)

  # Type II, test_intercepts = F
  LRTs <- compute_LRTs(fit$fit_obj,
                       ~ committee * emp_gender + (1 | id),
                       ~ committee * emp_gender + (1 | id),
                       dv = "assessment",
                       data = dat_exp_2,
                       type = 2,
                       mm)
  expect_equal(chisquares_two_factors_2[-c(1, 5)], as.numeric(LRTs$LRTs[, 4]), tolerance = 1e-4)

  expect_equal(stats::df.residual(LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 2)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)


  # Type III, test_intercepts = T
  LRTs <- compute_LRTs(fit$fit_obj,
                       ~ committee * emp_gender + (1 | id),
                       ~ committee * emp_gender + (1 | id),
                       dv = "assessment",
                       data = dat_exp_2,
                       type = 3,
                       mm,
                       test_intercepts = T)
  expect_equal(chisquares_two_factors_3, as.numeric(LRTs$LRTs[, 4]), tolerance = 1e-3)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[7]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[8]]), stats::df.residual(fit$fit_obj) + 1)

  # Type III, test_intercepts = T
  LRTs <- compute_LRTs(fit$fit_obj,
                       ~ committee * emp_gender + (1 | id),
                       ~ committee * emp_gender + (1 | id),
                       dv = "assessment",
                       data = dat_exp_2,
                       type = 3,
                       mm,
                       test_intercepts = F)
  expect_equal(chisquares_two_factors_3[-c(1, 5)], as.numeric(LRTs$LRTs[, 4]), tolerance = 1e-3)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[1]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[2]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[3]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[4]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[5]]), stats::df.residual(fit$fit_obj) + 1)
  expect_equal(stats::df.residual(LRTs$reduced_fits[[6]]), stats::df.residual(fit$fit_obj) + 1)
})



#------------------------------------------------------------------------------#
#### Factors with > 2 levels ####

test_that("compute_LRTs() works for factors with > 2 levels", {
  # Type II, test_intercepts = T
  form_mu <- ~ contingencies + (1 | id)
  form_lambda <- ~ contingencies + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mlsdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type_var = "status_ef",
                   data = dat_exp_2)

  mm <- construct_modelmatrices(form_mu, form_lambda,
                                dv = "assessment",
                                trial_type_var = "status",
                                data = dat_exp_2)

  # Type 2 - test_intercepts = T
  LRTs_2_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 2,
                                    mm,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_2_intercepts$LRTs[, 4]), chisquares_contingencies_2, tolerance = 1e-4)

  # Type 2 - test_intercepts = F
  LRTs_2 <- compute_LRTs(fit$fit_obj,
                         form_mu, form_lambda,
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 2,
                         mm,
                         test_intercepts = F)
  expect_equal(as.numeric(LRTs_2$LRTs[, 4]), chisquares_contingencies_2[c(2, 4)], tolerance = 1e-4)

  # Type 3 - test_intercepts = T
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = T)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chisquares_contingencies_3, tolerance = 1e-4)

  # Type 3 - test_intercepts = F
  LRTs_3_intercepts <- compute_LRTs(fit$fit_obj,
                                    form_mu, form_lambda,
                                    dv = "assessment",
                                    data = dat_exp_2,
                                    type = 3,
                                    mm,
                                    test_intercepts = F)
  expect_equal(as.numeric(LRTs_3_intercepts$LRTs[, 4]), chisquares_contingencies_3[c(2, 4)], tolerance = 1e-4)

})










