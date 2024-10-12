test_that("compute_tests() works with bootstraps on multiple cores", {
  # Type II, test_intercepts = T
  fit <- fit_mlsdt(formula_lambda = ~ committee * emp_gender + (1 | id),
                   formula_mu = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(4, "SOCK")
  boots <- compute_tests(fit$fit_obj,
                        ~ committee * emp_gender + (1 | id),
                        ~ committee * emp_gender + (1 | id),
                        dv = "assessment",
                        data = dat_exp_2,
                        type = 3,
                        tests = "bootstrap",
                        trial_type_var = "status_fac",
                        test_intercepts = T,
                        nsim = 4,
                        cl = cl)
  stopCluster(cl)
  expect_equal(boots$pb_objects[[1]]$parallel, T)
})

#------------------------------------------------------------------------------#
#### Parallelization LRTs ####

test_that("compute_tests() works with LRTs type 3 on multiple cores", {
  # Type II, test_intercepts = T
  fit <- fit_mlsdt(formula_lambda = ~ committee * emp_gender + (1 | id),
                   formula_mu = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(4, "SOCK")
  boots <- compute_tests(fit$fit_obj,
                         ~ committee * emp_gender + (1 | id),
                         ~ committee * emp_gender + (1 | id),
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 3,
                         tests = "bootstrap",
                         trial_type_var = "status_fac",
                         test_intercepts = T,
                         nsim = 4,
                         cl = cl)
  stopCluster(cl)

  expect_equal(boots$pb_objects[[1]]$parallel, T)
})


test_that("compute_tests() works with LRTs type 2 on multiple cores", {
  # Type II, test_intercepts = T
  fit <- fit_mlsdt(formula_lambda = ~ committee * emp_gender + (1 | id),
                   formula_mu = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(6, "SOCK")
  LRTs_par <- compute_tests(fit$fit_obj,
                         ~ committee * emp_gender + (1  | id),
                         ~ committee * emp_gender + (1 | id),
                         dv = "assessment",
                         data = dat_exp_2,
                         type = 2,
                         tests = "LRT",
                         trial_type_var = "status_fac",
                         test_intercepts = T,
                         cl = cl)
  stopCluster(cl)
  LRTs_seq <- compute_tests(fit$fit_obj,
                            ~ committee * emp_gender + (1 | id),
                            ~ committee * emp_gender + (1 | id),
                            dv = "assessment",
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            trial_type_var = "status_fac",
                            test_intercepts = T)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])
})
