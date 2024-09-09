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
                            (0 + mm[['rdm_lambda_ID']] + mm[['rdm_mu_ID']] | ID)"))
  )

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']]+
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
}
)


test_that("construct_glmer_formula() accepts different random effects grouping factors", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | ID),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv"
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] | VP) + (0 + mm[['rdm_mu_ID']] | ID)"))
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
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
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
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
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
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
})

test_that("construct_glmer_formula() makes a valid formula for uncorrelated random effects", {
  # Double bar notation
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 || ID),
      formula_lambda = ~ 1 + x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ 1 + x1 + (x1 || ID), ~1 + x1 + (x1 || ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_mu_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 2] | ID)"))
  )
  # No explicitly specified intercept
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (x1 || ID),
      formula_lambda = ~ x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x1 + (x1 | ID), ~x1 + (x1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_mu_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 2] | ID)"))
  )
  # Double bar notation + factor with multiple levels
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (x2 || ID),
      formula_lambda = ~ x2 + (x2 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (x2 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_lambda_ID']][, 3] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 1] | ID) + (0 + mm[['rdm_mu_ID']][, 2] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 3] | ID)"))
  )

  # Separate parentheses notation
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (1 | ID) + (x1 | ID),
      formula_lambda = ~ x1 + (1 | ID) + (x1 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_mu_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 2] | ID)"))
  )

  # Separate parentheses notation + multiple levels
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | ID) + (x2 | ID),
      formula_lambda = ~ x2 + (1 | ID) + (x2 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (x2 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_lambda_ID']][, 3] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 1] | ID) + (0 + mm[['rdm_mu_ID']][, 2] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 3] | ID)"))
  )

  # Different random-effects structures for lambda and mu
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | ID) + (x2 | ID),
      formula_lambda = ~ x2 + (1 | ID),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | ID), ~ x2 + (1 | ID), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 1] | ID) + (0 + mm[['rdm_mu_ID']][, 2] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 3] | ID)"))
  )

  # Different random factor
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x2 + (1 | rdm) + (x2 | rdm),
      formula_lambda = ~ x2 + (1 | rdm),
      dv = "dv",
      mm = construct_modelmatrices(~ x2 + (x2 | rdm), ~ x2 + (1 | rdm), data = internal_fake_data),
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_rdm']][, 1] | rdm) +
                            (0 + mm[['rdm_mu_rdm']][, 1] | rdm) + (0 + mm[['rdm_mu_rdm']][, 2] | rdm) +
                            (0 + mm[['rdm_mu_rdm']][, 3] | rdm)"))
  )
})

test_that("construct_glmer_formula() accepts different random effect identifiers.", {

  expect_equal(as.character(construct_glmer_formula(
    formula_mu = ~ x2 + (1 | rdm),
    formula_lambda = ~ x2 + (1 | rdm2),
    dv = "dv",
  )),
  as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_rdm2']] | rdm2) + (0 + mm[['rdm_mu_rdm']] | rdm)")))
})

test_that("construct_glmer_formula() removes model matrix if all columns are removed", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idc = Inf,
      remove_from_mu = T
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      param_idc = Inf,
      remove_from_mu = F
    )),
    as.character(as.formula("dv ~ 0 + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
})


#------------------------------------------------------------------------------#
#### correlate_sdt_params ####




test_that("construct_glmer_formula() handles correlate_sdt_params argument correctly", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      correlate_sdt_params = T,
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      correlate_sdt_params = F,
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] | VP) + (0 + mm[['rdm_mu_VP']] | VP)"))
  )

  # Uncorrelated terms
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 || ID),
      formula_lambda = ~ 1 + x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ 1 + x1 + (x1 | ID), ~ 1 + x1 + (x1 | ID), data = internal_fake_data),
      correlate_sdt_params = T
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_mu_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 2] | ID)"))
  )
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 || ID),
      formula_lambda = ~ 1 + x1 + (x1 || ID),
      dv = "dv",
      mm = construct_modelmatrices(~ 1 + x1 + (x1 | ID), ~ 1 + x1 + (x1 | ID), data = internal_fake_data),
      correlate_sdt_params = F
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_lambda_ID']][, 1] | ID) +
                            (0 + mm[['rdm_lambda_ID']][, 2] | ID) + (0 + mm[['rdm_mu_ID']][, 1] | ID) +
                            (0 + mm[['rdm_mu_ID']][, 2] | ID)"))
  )
})





#------------------------------------------------------------------------------#
#### Formulas with crossed random effects ####

test_that("construct_glmer_formula() handles crossed random effects properly", {

  expect_equal(
    # correlated mu and lambda
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (x1 | ID) + (x1 | stim),
      formula_lambda = ~ x2 + (x2 | ID) + (x2 | stim),
      dv = "y",
      correlate_sdt_params = F
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_ID']] | ID) + (0 + mm[['rdm_lambda_stim']] | stim) +
                            (0 + mm[['rdm_mu_ID']] | ID) + (0 + mm[['rdm_mu_stim']] | stim)")))


  # uncorrelated mu and lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + (x1 | ID) + (x1 | stim),
      formula_lambda = ~ x2 + (x2 | ID) + (x2 | stim),
      dv = "y",
      correlate_sdt_params = T
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_ID']] + mm[['rdm_mu_ID']] | ID) +
                          (0 + mm[['rdm_lambda_stim']] + mm[['rdm_mu_stim']] | stim)")))

  # uncorrelated random effects
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (committee || id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id) + (committee || stimulus),
      dv = "y",
      correlate_sdt_params = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (committee || id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee || id) + (committee || stimulus),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_lambda_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_lambda_stimulus']][, 2] | stimulus) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) + (0 + mm[['rdm_mu_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  # uncorrelated random effects that are different
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id) + (1 | stimulus),
      dv = "y",
      correlate_sdt_params = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee || id) + (1 | stimulus),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_lambda_stimulus']][, 1] | stimulus) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id) + (1 | stimulus),
      dv = "y",
      correlate_sdt_params = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee || id) + (1 | stimulus),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_lambda_stimulus']][, 1] | stimulus) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  # Different random-effects grouping factors for mu and lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']] | stimulus)")))
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))


  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = c(1, 2),
      remove_from_mu = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']][, -c(1, 2)] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = c(1, 3),
      remove_from_mu = F,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']][, -c(1, 3)] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))


}
)


#------------------------------------------------------------------------------#
#### Remove random effects ####

# TODO: handle uncorrelated random effects

test_that("construct_glmer_formula() removes indices from random effects", {
  # correlated random effects - remove from lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = c(1, 3),
      remove_from_mu = F,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, -c(1, 3)] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']] | stimulus)")))

  # correlated random effects - remove from mu
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = 2,
      remove_from_mu = T,
      remove_from_rdm = "stimulus",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, -2] | stimulus)")))

  # correlated random effects - don't correlated SDT params
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      correlate_sdt_params = F,
      remove_from_rdm = "VP",
      param_idc = 2
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']][, -2] | VP) + (0 + mm[['rdm_mu_VP']] | VP)"))
  )

  # correlated random effects - don't correlate SDT params
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP) + (x1 | ID),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      correlate_sdt_params = F,
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] | VP) + (0 + mm[['rdm_mu_VP']] | VP) + (0 + mm[['rdm_mu_ID']] | ID)"))
    )

  # remove whole model matrix when nothing is left in the model
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP) + (1 | ID),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      correlate_sdt_params = F,
      remove_from_rdm = "ID",
      remove_from_mu = T,
      param_idc = Inf
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] | VP) + (0 + mm[['rdm_mu_VP']] | VP)"))
  )
})


# TODO: uncorrelated random effects
test_that("construct_glmer_formula() removes indices from random effects", {
  # correlated random effects - remove from lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = c(1, 3),
      remove_from_mu = F,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, -c(1, 3)] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']] | stimulus)")))

  # correlated random effects - remove from mu
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = c(1, 3),
      remove_from_mu = T,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] + mm[['rdm_mu_id']][, -c(1, 3)] | id) +
                            (0 + mm[['rdm_mu_stimulus']] | stimulus)")))

  # correlated random effects - remove from mu (stimulus)
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = c(1, 3),
      remove_from_mu = T,
      remove_from_rdm = "stimulus",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, -c(1, 3)] | stimulus)")))

  # correlated random effects - remove whole matrix when nothing is left
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      param_idc = Inf,
      remove_from_mu = T,
      remove_from_rdm = "stimulus",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] + mm[['rdm_mu_id']] | id)")))

  # correlated random effects - uncorrelated mu and lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = c(1, 3),
      remove_from_mu = T,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] | id) + (0 + mm[['rdm_mu_id']][, -c(1, 3)] | id) +
                            (0 + mm[['rdm_mu_stimulus']] | stimulus)")))


  # correlated random effects - remove whole matrix when nothing is left
  # uncorrelated mu and lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = Inf,
      remove_from_mu = T,
      remove_from_rdm = "stimulus",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']] | id) +(0 + mm[['rdm_mu_id']] | id)")))
})





test_that("construct_glmer_formula() removes indices from random effects (uncorrelated)", {
  # uncorrelated random effects - remove from lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id),
      dv = "y",
      param_idc = 1,
      remove_from_mu = F,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee || id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  # uncorrelated random effects - remove from mu
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id),
      dv = "y",
      param_idc = 1,
      remove_from_mu = T,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee || stimulus),
                                   formula_lambda = ~ committee + (committee || id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  # uncorrelated random effects - remove whole matrix when nothing is left
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = Inf,
      remove_from_mu = T,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  # uncorrelated random effects - remove whole matrix when nothing is left (lambda)
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id),
      dv = "y",
      correlate_sdt_params = F,
      param_idc = Inf,
      remove_from_mu = F,
      remove_from_rdm = "id",
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))
})


#------------------------------------------------------------------------------#
#### No random effects ####

test_that("construct_glmer_formula() makes the correct formula without random effects", {
  # -> We don't need mm here
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ x1 + x2,
      formula_lambda = ~ x1 + x2,
      dv = "y"
    )),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']]"))
  )
}
)

#------------------------------------------------------------------------------#
#### remove_correlations argument ####

test_that("construct_glmer_formula() removes correlations from correlated formula", {
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      remove_correlations = T,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      remove_correlations = T,
      remove_from_mu = F,
      remove_from_rdm = "id",
      param_idc = 1,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      remove_correlations = T,
      remove_from_mu = T,
      remove_from_rdm = "stimulus",
      param_idc = Inf,
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id)")))
}
)
