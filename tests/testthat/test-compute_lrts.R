
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

test_that("compute_LRTs() throws a message when trying to test an intercept without predictors", {
  fit <- fit_mlsdt(~ 1 + (x1 || ID), ~ 1 + (x1 || ID), dv = "y", data = internal_sdt_data)$fit_obj
  mm <- construct_modelmatrices(~ 1 + (x1 || ID), ~ 1 + (x1 || ID), data = internal_sdt_data)
  expect_message(compute_LRTs(fit, ~ 1 + (x1 || ID), ~ 1 + (x1 || ID), dv = "y", data = internal_sdt_data,
                              mm  = mm, test_intercepts = T))
  #expect_equal(compute_LRTs(fit, ~ 1 + (x1 || ID), ~ 0 + (x1 || ID), dv = "y", data = internal_sdt_data,
  #                          mm  = mm, test_intercepts = T),
  #             NULL)

})


#------------------------------------------------------------------------------#
#### Type II SS ####

test_that("compute_LRTs() Type II works with a standard two-factorial design without test_intercepts", {
  # Same for the uncorrelated model
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
  expect_equal(chisquares_two_factors, unlist(LRTs$LRTs[, 4]), tolerance = 1e-4)
})




test_that("compute_LRTs() Type II works with a standard two-factorial design with test_intercepts", {
  # Same for the uncorrelated model
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
                       mm)
  expect_equal(chisquares_two_factors[-c(1, 5)], unlist(LRTs$LRTs[, 4]), tolerance = 1e-4)
})


# Appropriate error messages










