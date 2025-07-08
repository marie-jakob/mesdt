#------------------------------------------------------------------------------#
#### Removal from fixed effects ####

test_that("construct_glmer_formula() makes a valid reduced formula", {

  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      to_remove = list("mu" = 2)
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']][, -2] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ 1 + x1 + (x1 | VP),
      formula_lambda = ~ 1 + x2 + (x2 | VP),
      dv = "dv",
      to_remove = list("lambda" = 3)
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
      to_remove = list("mu" = which(c(1, 3, 1) == 1))
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']][, -c(1, 3)] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)"))
  )
})


#------------------------------------------------------------------------------#
#### Removal from random effects ####


test_that("construct_glmer_formula() removes indices from random effects", {
  # correlated random effects - remove from lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list("rdm_lambda_id" = c(1, 3)),
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
      to_remove = list("rdm_mu_stimulus" = c(2)),
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
      to_remove = list("rdm_lambda_VP" = 2)
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
      to_remove = list("rdm_mu_ID" = Inf)
    )),
    as.character(as.formula("dv ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] | VP) + (0 + mm[['rdm_mu_VP']] | VP)"))
  )
})


test_that("construct_glmer_formula() removes indices from random effects", {
  # correlated random effects - remove from lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list("rdm_lambda_id" = c(1, 3)),
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
      to_remove = list("rdm_mu_id" = c(1, 3)),
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
      to_remove = list("rdm_mu_stimulus" = c(1, 3)),
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
      to_remove = list("rdm_mu_stimulus" = Inf),
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
      to_remove = list("rdm_mu_id" = c(1, 3)),
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
      to_remove = list("rdm_mu_stimulus" = Inf),
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
      to_remove = list("rdm_lambda_id" = 1),
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
      to_remove = list("rdm_mu_id" = 1),
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
      to_remove = list("rdm_mu_id" = Inf),
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
      to_remove = list("rdm_lambda_id" = Inf),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus) + (0 + mm[['rdm_mu_stimulus']][, 2] | stimulus)")))
})


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
      to_remove = list("rdm_lambda_id" = 1),
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
      to_remove = list("rdm_mu_stimulus" = Inf),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_id']][, 1] | id) + (0 + mm[['rdm_lambda_id']][, 2] | id) +
                            (0 + mm[['rdm_mu_id']][, 1] | id)")))
}
)

#------------------------------------------------------------------------------#
#### Removal from multiple parts ####

test_that("construct_glmer_formula() removes indices from multiple parts at once", {

  # remove a mix of everything - correlated
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list(
        "lambda" = 1,
        "mu" = 2,
        "rdm_lambda_id" = 1,
        "rdm_mu_stimulus" = 2),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']][, -1] + mm[['mu']][, -2] +
                            (0 + mm[['rdm_lambda_id']][, -1] + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, -2] | stimulus)")))

  # remove whole matrices
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list(
        "rdm_lambda_id" = Inf,
        "rdm_mu_stimulus" = 1),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']] + (0 + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, -1] | stimulus)")))

  # edge cases: all random effects removed -> should be possible
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list(
        "rdm_lambda_id" = Inf,
        "rdm_mu_stimulus" = Inf,
        "rdm_mu_id" = Inf),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']] + mm[['mu']]")))

  # all fixed effects removed -> doesn't make sense, throw a warning
  expect_message(
    construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list(
        "lambda" = Inf,
        "mu" = Inf),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")))
  expect_equal(
    construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = T,
      to_remove = list(
        "lambda" = Inf,
        "mu" = Inf),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac")), NULL)

  # remove a mix of everything - uncorrelated mu and lambda
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee | stimulus),
      formula_lambda = ~ committee + (committee | id),
      dv = "y",
      correlate_sdt_params = F,
      to_remove = list(
        "lambda" = 1,
        "mu" = 2,
        "rdm_lambda_id" = 1,
        "rdm_mu_stimulus" = 2),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']][, -1] + mm[['mu']][, -2] +
                            (0 + mm[['rdm_lambda_id']][, -1] | id) + (0 + mm[['rdm_mu_id']] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, -2] | stimulus)")))

  # remove a mix of everything - everything uncorrelated
  expect_equal(
    as.character(construct_glmer_formula(
      formula_mu = ~ committee + (1 | id) + (committee || stimulus),
      formula_lambda = ~ committee + (committee || id),
      dv = "y",
      correlate_sdt_params = F,
      to_remove = list(
        "lambda" = 1,
        "mu" = 2,
        "rdm_lambda_id" = 1,
        "rdm_mu_stimulus" = 2),
      mm = construct_modelmatrices(formula_mu = ~ committee + (1 | id) + (committee | stimulus),
                                   formula_lambda = ~ committee + (committee | id),
                                   data = dat_exp_2,
                                   trial_type_var = "status_fac"))),
    as.character(as.formula("y ~ 0 + mm[['lambda']][, -1] + mm[['mu']][, -2] +
                            (0 + mm[['rdm_lambda_id']][, 2] | id) + (0 + mm[['rdm_mu_id']][, 1] | id) +
                            (0 + mm[['rdm_mu_stimulus']][, 1] | stimulus)")))

})


