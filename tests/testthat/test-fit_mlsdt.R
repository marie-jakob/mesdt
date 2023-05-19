#------------------------------------------------------------------------------#
#### construct_glmer_formula() ####


test_that("construct_glmer_formula() makes the correct formula", {
  expect_equal(
    as.character(construct_glmer_formula(
    formula_mu = ~ x1 + (1 | ID),
    formula_lambda = ~ x2 + (1 | ID),
    dv = "y"
    )),
    as.character(as.formula("y ~ 0 + modeldata_lambda + modeldata_mu + (0 + modeldata_random_lambda + modeldata_random_mu | ID)"))
  )

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    )),
    as.character(as.formula("dv ~ 0 + modeldata_lambda + modeldata_mu + (0 + modeldata_random_lambda + modeldata_random_mu | VP)"))
  )

  expect_equal(
    construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | ID),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    ),
    NULL
  )
}
)

#------------------------------------------------------------------------------#
#### construct_modeldata() ####



test_that("construct_modeldata() constructs valid modeldata", {
  expect_equal(
    construct_modeldata(),
    as.character(as.formula("y ~ 0 + modeldata_lambda + modeldata_mu + (0 + modeldata_random_lambda + modeldata_random_mu | ID)"))
  )
}


test_that("construct_modeldata() works for mu", {
  expect_equal(
    construct_modeldata()),
    as.character(as.formula("y ~ 0 + modeldata_lambda + modeldata_mu + (0 + modeldata_random_lambda + modeldata_random_mu | ID)"))
  )
}

