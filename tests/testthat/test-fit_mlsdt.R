#------------------------------------------------------------------------------#
#### construct_glmer_formula() ####

test_that("construct_glmer_formula() makes the correct formula", {
  # -> We don't need mm here
  expect_equal(
    as.character(construct_glmer_formula(
    formula_mu = ~ x1 + (1 | ID),
    formula_lambda = ~ x2 + (1 | ID),
    dv = "y"
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']]+
                            (0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | ID)"))
  )

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']]+
                            (0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | VP)"))
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
      param_idc = 2,
      remove_from_mu = T
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']][, -2] +
                            (0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | VP)"))
  )
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idc = 3,
      remove_from_mu = F
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']][, -3] + mm[['mu']] +
                            (0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | VP)"))
  )
})

test_that("construct_glmer_formula() makes a valid reduced formula for a vector of indices", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idc = which(c(1, 3, 1) == 1),
      remove_from_mu = T
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']][, -c(1, 3)] +
                            (0 + mm[['rdm_lambda']] + mm[['rdm_mu']] | VP)"))
  )
})

test_that("construct_glmer_formula() makes a valid formula for uncorrelated random effects", {
  # Double bar notation
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 || ID),
      formula_lambda = ~ 1 + x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ 1 + x1 + (x1 | ID), ~1 + x1 + (x1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_lambda']][, 2] | ID) + (0 + mm[['rdm_mu']][, 1] | ID) +
                            (0 + mm[['rdm_mu']][, 2] | ID)"))
  )
  # No explicitly specified intercept
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (x1 || ID),
      formula_lambda = ~ x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x1 + (x1 | ID), ~x1 + (x1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_lambda']][, 2] | ID) + (0 + mm[['rdm_mu']][, 1] | ID) +
                            (0 + mm[['rdm_mu']][, 2] | ID)"))
  )
  # Double bar notation + factor with multiple levels
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (x2 || ID),
      formula_lambda = ~ x2 + (x2 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (x2 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_lambda']][, 2] | ID) + (0 + mm[['rdm_lambda']][, 3] | ID) +
                            (0 + mm[['rdm_mu']][, 1] | ID) + (0 + mm[['rdm_mu']][, 2] | ID) +
                            (0 + mm[['rdm_mu']][, 3] | ID)"))
  )

  # Separate parentheses notation
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (1 | ID) + (x1 | ID),
      formula_lambda = ~ x1 + (1 | ID) + (x1 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_lambda']][, 2] | ID) + (0 + mm[['rdm_mu']][, 1] | ID) +
                            (0 + mm[['rdm_mu']][, 2] | ID)"))
  )

  # Separate parentheses notation + multiple levels
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | ID) + (x2 | ID),
      formula_lambda = ~ x2 + (1 | ID) + (x2 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (x2 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_lambda']][, 2] | ID) + (0 + mm[['rdm_lambda']][, 3] | ID) +
                            (0 + mm[['rdm_mu']][, 1] | ID) + (0 + mm[['rdm_mu']][, 2] | ID) +
                            (0 + mm[['rdm_mu']][, 3] | ID)"))
  )

  # Different random-effects structures for lambda and mu
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | ID) + (x2 | ID),
      formula_lambda = ~ x2 + (1 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | ID) +
                            (0 + mm[['rdm_mu']][, 1] | ID) + (0 + mm[['rdm_mu']][, 2] | ID) +
                            (0 + mm[['rdm_mu']][, 3] | ID)"))
  )

  # Different random factor
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | rdm) + (x2 | rdm),
      formula_lambda = ~ x2 + (1 | rdm),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | rdm), ~ x2 + (1 | rdm), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda']][, 1] | rdm) +
                            (0 + mm[['rdm_mu']][, 1] | rdm) + (0 + mm[['rdm_mu']][, 2] | rdm) +
                            (0 + mm[['rdm_mu']][, 3] | rdm)"))
  )
})


test_that("construct_glmer_formula() does not accept different random effect identifiers.", {
  expect_message(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2 + (1 | rdm2),
    dv = "dv",
  ))

  expect_equal(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2 + (1 | rdm2),
    dv = "dv",
  ),
  NULL)
})

test_that("construct_glmer_formula() does not accept formulas without random effects.", {
  expect_message(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2,
    dv = "dv",
  ))

  expect_equal(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2,
    dv = "dv",
  ), NULL)

  expect_message(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2,
    dv = "dv",
  ))
  expect_equal(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2,
    dv = "dv",
  ), NULL)
})

#------------------------------------------------------------------------------#
#### construct_modelmatrices() ####

test_that("construct_modelmatrices() constructs valid mm for a single predictor and a random intercept", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 + (1 | ID),
                        formula_lambda = ~ x1 + (1 | ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1, data = internal_fake_data),
         "rdm_mu" = stats::model.matrix(~ 1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda" = stats::model.matrix(~ 1, data = internal_fake_data)
  ))
})

test_that("construct_modelmatrices() constructs valid mm for multiple predictors and random slopes", {
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (x1 | ID),
                        formula_lambda = ~ x1 + x2 + (x2 | ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda" = stats::model.matrix(~ x2, data = internal_fake_data)
  ))
}
)

test_that("construct_modelmatrices() works for suppressed correlations", {
  # -> mm should not be affected, only formula

  # "||" syntax
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (x1 || ID),
                        formula_lambda = ~ x1 + x2 + (x2 || ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda" = stats::model.matrix(~ x2, data = internal_fake_data)
    ))

  # separated random effects syntax
  expect_equal(
    construct_modelmatrices(formula_mu = ~ x1 * x2 + (1 | ID) + (0 + x1 | ID),
                        formula_lambda = ~ x1 + x2 + (1 | ID) + (0 + x2 | ID),
                        dv = "y",
                        data = internal_fake_data),
    list("mu" = stats::model.matrix(~ x1 * x2, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "lambda" = stats::model.matrix(~ x1 + x2, data = internal_fake_data),
         "rdm_mu" = stats::model.matrix(~ x1, data = internal_fake_data) *
           stats::model.matrix(~ trial_type, data = internal_fake_data)[, 2] * 0.5,
         "rdm_lambda" = stats::model.matrix(~ x2, data = internal_fake_data)
    ))
})

#------------------------------------------------------------------------------#
#### fit_mlsdt() ####

# use a saved model for this

test_that("fit_mlsdt() estimates the correct model", {
  fit <- fit_mlsdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data)$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(model_test)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(model_test)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test))))

  # fixed effects estimates
  expect_equal(unname(fixef(fit))[1:2], unname(fixef(model_test))[1:2], tolerance = 1e-4)
  # mu fixef effects
  expect_equal(unname(fixef(fit))[3:4], unname(fixef(model_test))[3:4] * 2, tolerance = 1e-4)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1:2, 1:2], unname(vcov(model_test))[1:2, 1:2], tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[3:4, 3:4], unname(vcov(model_test))[3:4, 3:4] * 4, tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[1:2, 3:4], unname(vcov(model_test))[1:2, 3:4] * 2, tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(as.data.frame(VarCorr(fit))$vcov[1:2], as.data.frame(VarCorr(model_test))$vcov[1:2], tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[3:4], as.data.frame(VarCorr(model_test))$vcov[3:4] * 4, tolerance = 1e-3)

  # random effects correlations
  expect_equal(as.data.frame(VarCorr(fit))$sdcor[5:10],
               as.data.frame(VarCorr(model_test))$sdcor[5:10],
               tolerance = 1e-3)
}
)


test_that("fit_mlsdt() works for uncorrelated random effects (|| notation)", {
  fit <- fit_mlsdt(~ 1 + x1 + (1 + x1 || ID), ~ 1 + x1 + (1 + x1 || ID), dv = "y", data = internal_sdt_data)$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(model_test_uncor)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(model_test_uncor)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test_uncor))))

  # fixed effects estimates
  expect_equal(unname(fixef(fit))[1:2], unname(fixef(model_test_uncor))[1:2], tolerance = 1e-4)
  # mu fixef effects
  expect_equal(unname(fixef(fit))[3:4], unname(fixef(model_test_uncor))[3:4] * 2, tolerance = 1e-4)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1:2, 1:2], unname(vcov(model_test_uncor))[1:2, 1:2], tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[3:4, 3:4], unname(vcov(model_test_uncor))[3:4, 3:4] * 4, tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[1:2, 3:4], unname(vcov(model_test_uncor))[1:2, 3:4] * 2, tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(as.data.frame(VarCorr(fit))$vcov[1:2], as.data.frame(VarCorr(model_test_uncor))$vcov[1:2], tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[3:4], as.data.frame(VarCorr(model_test_uncor))$vcov[3:4] * 4, tolerance = 1e-3)

  # random effects correlations
  expect_equal(as.data.frame(VarCorr(fit))$sdcor[5:10],
               as.data.frame(VarCorr(model_test_uncor))$sdcor[5:10],
               tolerance = 1e-3)
}
)

test_that("fit_mlsdt() notifies the user that only uncorrelated or correlated
          random effects are possible at the moment, iff specified.", {
  expect_message(fit_mlsdt(~ 1 + x1 + (1 + x1 || ID), ~ 1 + x1 + (1 + x1 | ID), dv = "y", data = internal_sdt_data))
  expect_message(fit_mlsdt(~ 1 + x1 + (1 | ID) + (x1 | ID), ~ 1 + x1 + (1 + x1 || ID), dv = "y", data = internal_sdt_data))

  expect_no_message(fit_mlsdt(~ 1 + x1 + (x1 | ID), ~ 1 + x1 + (x1 | ID), dv = "y", data = internal_sdt_data))
  expect_no_message(fit_mlsdt(~ 1 + x1 + (1 | ID), ~ 1 + x1 + (1 | ID), dv = "y", data = internal_sdt_data))

})

