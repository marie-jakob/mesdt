#------------------------------------------------------------------------------#
#### construct_glmer_formula() ####

test_that("construct_glmer_formula() makes the correct formula", {
  expect_equal(
    as.character(construct_glmer_formula(
    formula_mu = ~ x1 + (1 | ID),
    formula_lambda = ~ x2 + (1 | ID),
    dv = "y"
    )),
    as.character(as.formula("y ~ 0 + modeldata[['lambda']] + modeldata[['mu']]+
                            (0 + modeldata[['random_lambda']] + modeldata[['random_mu']] | ID)"))
  )

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    )),
    as.character(as.formula("dv ~ 0 + modeldata[['lambda']] + modeldata[['mu']]+
                            (0 + modeldata[['random_lambda']] + modeldata[['random_mu']] | VP)"))
  )
}
)

test_that("construct_glmer_formula() handles wrong input properly", {
  expect_message(construct_glmer_formula(
    formula_mu = ~ 1 + x1 + (x1 | ID),
    formula_lambda = ~ 1 + x2 + (x2 | VP),
    dv = "dv"
  ))
  expect_equal(
    construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | ID),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    ),
    NULL
  )
})

test_that("construct_glmer_formula() makes a valid reduced formula", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idx = 2,
      remove_from_mu = T
    )),
    as.character(as.formula("dv ~ 0 + modeldata[['lambda']] + modeldata[['mu']][, -2] +
                            (0 + modeldata[['random_lambda']] + modeldata[['random_mu']] | VP)"))
  )
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idx = 3,
      remove_from_mu = F
    )),
    as.character(as.formula("dv ~ 0 + modeldata[['lambda']][, -3] + modeldata[['mu']] +
                            (0 + modeldata[['random_lambda']] + modeldata[['random_mu']] | VP)"))
  )
})

#------------------------------------------------------------------------------#
#### construct_modeldata() ####

test_that("construct_modeldata() constructs valid modeldata for a single predictor and a random intercept", {
  expect_equal(
    construct_modeldata(formula_mu = ~ x1 + (1 | ID),
                        formula_lambda = ~ x1 + (1 | ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2],
         "lambda" = stats::model.matrix(~ x1, data = internal_fake_data),
         "random_mu" = stats::model.matrix(~ 1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2],
         "random_lambda" = stats::model.matrix(~ 1, data = internal_fake_data)
  ))
})

test_that("construct_modeldata() constructs valid modeldata for multiple predictors and random slopes", {
  expect_equal(
    construct_modeldata(formula_mu = ~ x1 * x2 + (x1 | ID),
                        formula_lambda = ~ x1 + x2 + (x2 | ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2],
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "random_mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2],
         "random_lambda" = stats::model.matrix(~ x2, data = internal_fake_data)
  ))
}
)
