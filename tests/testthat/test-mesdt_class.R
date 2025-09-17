
test_that("simulate() does something sensible", {
  skip_if_not_installed("glmmTMB")
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                   trial_type = "status_fac")

  test_sim <- simulate(fit, nsim = 10, seed = 3)
  expect_equal(ncol(test_sim), 10)

  # check if seed is used correctly
  test_sim_1 <- simulate(fit, nsim = 1, seed = 3)
  test_sim_2 <- simulate(fit, nsim = 1, seed = 3)
  expect_equal(test_sim_1, test_sim_2)

  options("mesdt.backend" = "glmmTMB")
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                   trial_type = "status_fac")

  test_sim <- simulate(fit, nsim = 10, seed = 3)
  expect_equal(ncol(test_sim), 10)
  # check if seed is used correctly
  test_sim_1 <- simulate(fit, nsim = 1, seed = 3)
  test_sim_2 <- simulate(fit, nsim = 1, seed = 3)
  expect_equal(test_sim_1, test_sim_2)


  options("mesdt.backend" = "lme4")
})




test_that("summary() works for factors with > 2 levels", {
  form_mu <- ~ contingencies + (1 | id)
  form_lambda <- ~ contingencies + (1 | id)
  # Same for the uncorrelated model
  fit <- fit_mesdt(form_mu, form_lambda,
                   dv = "assessment",
                   trial_type = "status_fac",
                   data = dat_exp_2)
  s <- summary(fit)


  expect_equal(unname(s$c_coef[, 1]), unname(fixef(fit$fit_obj)[1:3] * -1))
  expect_equal(unname(s$d_coef[, 1]), unname(fixef(fit$fit_obj)[4:6]))

})

