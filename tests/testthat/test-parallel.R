

options("mesdt.backend" = "lme4")

test_that("compute_tests() works with bootstraps on multiple cores", {
  library(parallel)
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(4, "SOCK")
  boots_cl <- compute_tests(fit,
                        data = dat_exp_2,
                        type = 3,
                        tests = "bootstrap",
                        test_intercepts = T,
                        nsim = 4,
                        cl = cl)
  parallel::stopCluster(cl)
  expect_equal(boots_cl$pb_objects[[1]]$parallel, T)
})

#------------------------------------------------------------------------------#
#### Parallelization LRTs ####

test_that("compute_tests() works with LRTs type 3 on multiple cores", {
  library(parallel)
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(6, "SOCK")
  parallel::clusterEvalQ(cl = cl, {options("mesdt.backend" = "lme4")})
  LRTs_par <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 3,
                            tests = "LRT",
                            test_intercepts = T,
                            cl = cl)
  parallel::stopCluster(cl)
  options("mesdt.backend" = "lme4")
  LRTs_seq <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 3,
                            tests = "LRT",
                            test_intercepts = T)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])
})


test_that("compute_tests() works with LRTs type 2 on multiple cores", {
  library(parallel)
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(6, "SOCK")
  # parallel::clusterEvalQ(cl = cl, {options("mesdt.backend" = "lme4")})
  LRTs_par <- compute_tests(fit,
                         data = dat_exp_2,
                         type = 2,
                         tests = "LRT",
                         test_intercepts = T,
                         cl = cl)
  parallel::stopCluster(cl)
  options("mesdt.backend" = "lme4")
  LRTs_seq <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            test_intercepts = T)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])
})

test_that("compute_tests() works with LRTs using glmmTMB as backend", {
  library(parallel)
  options("mesdt.backend" = "glmmTMB")
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  cl <- parallel::makeCluster(6, "SOCK")

  LRTs_par <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            test_intercepts = T,
                            cl = cl)
  parallel::stopCluster(cl)

  LRTs_seq <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            test_intercepts = T)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])

  cl <- parallel::makeCluster(6, "SOCK")

  LRTs_par <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 3,
                            tests = "LRT",
                            test_intercepts = T,
                            cl = cl)
  parallel::stopCluster(cl)

  LRTs_seq <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 3,
                            tests = "LRT",
                            test_intercepts = T)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])
})


test_that("compute_tests() sets the correct backend", {
  library(parallel)
  options("mesdt.backend" = "glmmTMB")
  # Type II, test_intercepts = T
  fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                   discriminability = ~ committee * emp_gender + (1 | id),
                   dv = "assessment",
                   trial_type_var = "status_fac",
                   data = dat_exp_2)
  LRTs_seq <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            test_intercepts = T)
  options("mesdt.backend" = "lme4")
  cl <- parallel::makeCluster(6, "SOCK")
  LRTs_par <- compute_tests(fit,
                            data = dat_exp_2,
                            type = 2,
                            tests = "LRT",
                            test_intercepts = T,
                            cl = cl)
  parallel::stopCluster(cl)
  expect_equal(LRTs_par$LRTs[, 4], LRTs_seq$LRTs[, 4])

})
