
#------------------------------------------------------------------------------#
#### construct_modelmatrices() ####

test_that("construct_modelmatrices() constructs valid mm for a single predictor and a random intercept", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 + (1 | ID),
                            formula_lambda = ~ x1 + (1 | ID),
                            data = internal_fake_data)[["mm"]],
    list("mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1, data = internal_fake_data),
         "rdm_mu_ID" = stats::model.matrix(~ 1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda_ID" = stats::model.matrix(~ 1, data = internal_fake_data)
    ))
})

test_that("construct_modelmatrices() constructs valid mm for multiple predictors and random slopes", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (x1 | ID),
                            formula_lambda = ~ x1 + x2 + (x2 | ID),
                            data = internal_fake_data)[["mm"]],
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu_ID" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda_ID" = stats::model.matrix(~ x2, data = internal_fake_data)
    ))
}
)

test_that("construct_modelmatrices() works for suppressed correlations", {
  # -> mm should not be affected, only formula
  # "||" syntax
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (x1 || ID),
                            formula_lambda = ~ x1 + x2 + (x2 || ID),
                            data = internal_fake_data)[["mm"]],
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu_ID" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda_ID" = stats::model.matrix(~ x2, data = internal_fake_data)
    ))

  # separated random effects syntax
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (1 | ID) + (0 + x1 | ID),
                            formula_lambda = ~ x1 + x2 + (1 | ID) + (0 + x2 | ID),
                            data = internal_fake_data)[["mm"]],
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu_ID" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda_ID" = stats::model.matrix(~ x2, data = internal_fake_data)
    ))
})

#------------------------------------------------------------------------------#
#### Crossed random effects ####


test_that("construct_modelmatrices() constructs valid mm for crossed random effects", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ 1 + (1 | id) + (1 | file_name),
                            formula_lambda = ~ 1 + (1 | id) + (1 | file_name),
                            trial_type_var = "status_fac",
                            data = dat_exp_2)[["mm"]],
    list("mu" = stats::model.matrix(~ 1, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ 1, data = dat_exp_2),
         "rdm_mu_id" = stats::model.matrix(~ 1, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_mu_file_name" = stats::model.matrix(~ 1, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_lambda_id" = stats::model.matrix(~ 1, data = dat_exp_2),
         "rdm_lambda_file_name" = stats::model.matrix(~ 1, data = dat_exp_2)
    ))
})


test_that("construct_modelmatrices() constructs valid mm for crossed random effects with one predictor and random slopes", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ committee + (committee | id) + (committee | file_name),
                            formula_lambda = ~ committee + (committee | id) + (committee | file_name),
                            trial_type_var = "status_fac",
                            data = dat_exp_2)[["mm"]],
    list("mu" = stats::model.matrix(~ committee, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ committee, data = dat_exp_2),
         "rdm_mu_id" = stats::model.matrix(~ committee, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_mu_file_name" = stats::model.matrix(~ committee, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_lambda_id" = stats::model.matrix(~ committee, data = dat_exp_2),
         "rdm_lambda_file_name" = stats::model.matrix(~ committee, data = dat_exp_2)
    ))
})

test_that("construct_modelmatrices() constructs valid mm for different random effects for mu and lambda", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ 1 + (emp_gender | id),
                            formula_lambda = ~ 1 + (committee | file_name),
                            trial_type_var = "status_fac",
                            data = dat_exp_2)[["mm"]],
    list("mu" = stats::model.matrix(~ 1, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ 1, data = dat_exp_2),
         "rdm_mu_id" = stats::model.matrix(~ emp_gender, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_lambda_file_name" = stats::model.matrix(~ committee, data = dat_exp_2)
    ))
})



test_that("construct_modelmatrices() constructs valid mm for removed intercept from fixed effects", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ 0 + emp_gender + (emp_gender | id),
                            formula_lambda = ~ 0 + emp_gender + (committee | file_name),
                            trial_type_var = "status_fac",
                            data = dat_exp_2)[["mm"]],
    list("mu" = stats::model.matrix(~ 0 + emp_gender, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ 0 + emp_gender, data = dat_exp_2),
         "rdm_mu_id" = stats::model.matrix(~ emp_gender, data = dat_exp_2) *
           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
         "rdm_lambda_file_name" = stats::model.matrix(~ committee, data = dat_exp_2)
    ))
})



# TODO: fix that
#test_that("construct_modelmatrices() constructs valid mm for removed intercept from random effects", {
#  expect_equal(
#    construct_modelmatrices(formula_mu = ~ 1 + (0 + emp_gender | id),
#                            formula_lambda = ~ 1 + (0 + committee | file_name),
#                            dv = "y",
#                            trial_type_var = "status_fac",
#                            data = dat_exp_2),
#    list("mu" = stats::model.matrix(~ 1, data = dat_exp_2) *
#           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
#         "lambda" = stats::model.matrix(~ 1, data = dat_exp_2),
#         "rdm_mu_id" = stats::model.matrix(~ 0 + emp_gender, data = dat_exp_2) *
#           stats::model.matrix(~ status_fac, data = dat_exp_2)[, 2] * 0.5,
#         "rdm_lambda_file_name" = stats::model.matrix(~ 0 + committee, data = dat_exp_2)
#    ))
#})


