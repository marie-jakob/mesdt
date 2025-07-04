#------------------------------------------------------------------------------#
#### fit_mesdt() ####


for (backend in c("glmmTMB", "lme4")) {
  options("mesdt.backend" = backend)


  test_that("compute_tests() compares the correct models for bootstrap tests of fixed effects", {

    fit <- fit_mesdt(bias = ~ committee * emp_gender + (1 | id),
                     discriminability = ~ committee * emp_gender_ef + (1 | id),
                     dv = "assessment",
                     trial_type_var = "status_fac",
                     data = dat_exp_2)
    for (type_tmp in c(2, 3)) {
      for (inter_tmp in c(T, F)) {
        print(inter_tmp)
        print(type_tmp)
        suppressWarnings(boot <- compute_tests(fit,
                              tests = "bootstrap",
                              type = type_tmp,
                              nsim = 1,
                              test_intercepts = inter_tmp))
        expect_equal(boot$pb_test_results[, 1], boot$LRT_results[, 4])
      }
    }
  })

  test_that("compute_tests() compares the correct models for bootstrap tests of factors with multiple levels", {

    form_mu <- ~ contingencies + (1 | id)
    form_lambda <- ~ contingencies + (1 | id)
    # Same for the uncorrelated model
    fit <- fit_mesdt(form_mu, form_lambda,
                     dv = "assessment",
                     trial_type_var = "status_fac",
                     data = dat_exp_2)
    for (inter_tmp in c(T, F)) {
      for (type_tmp in c(2, 3)) {
        print(inter_tmp)
        print(type_tmp)
        suppressWarnings(boot <- compute_tests(fit,
                              test_intercepts = inter_tmp,
                              test_ran_ef = F,
                              type = type_tmp,
                              tests = "bootstrap",
                              nsim = 1))
        expect_equal(boot$pb_test_results[, 1], boot$LRT_results[, 4])
      }
    }
  })


  test_that("compute_tests() compares the correct models for bootstrap tests of random effects", {

    fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                     trial_type_var = "status_fac")
    for (inter_tmp in c(T, F)) {
      print(inter_tmp)
      suppressWarnings(boot <- compute_tests(fit,
                           test_intercepts = T,
                           test_ran_ef = inter_tmp,
                           type = 3,
                           tests = "bootstrap",
                           nsim = 1))
      expect_equal(boot$pb_test_results[, 1], boot$LRT_results[, 4])
    }
  })

  test_that("compute_tests() uses a given cluster for bootstrap tests", {

    skip_if_not_installed("parallel")
    library(parallel)

    fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                     trial_type_var = "status_fac")
    cl <- makeCluster(8, type = "SOCK")
    suppressWarnings(boot <- compute_tests(fit,
                                           test_intercepts = T,
                                           test_ran_ef = T,
                                           type = 3,
                                           tests = "bootstrap",
                                           nsim = 8,
                                           cl = cl))

      expect_equal(unname(boot$pb_objects[[1]]$parallel), T)
      stopCluster(cl)
  })
}


options("mesdt.backend" = "lme4")


#------------------------------------------------------------------------------#
#### Section ####

test_that("compute_tests() uses the correct seed for bootstrapping", {

  skip_if_not_installed("parallel")
  library(parallel)

  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name), ~ committee + (committee | id), dv = "assessment", data = dat_exp_2,
                   trial_type_var = "status_fac")
  cl <- makeCluster(8, type = "SOCK")
  suppressWarnings(boot_1 <- compute_tests(fit,
                                         test_intercepts = T,
                                         test_ran_ef = T,
                                         type = 3,
                                         tests = "bootstrap",
                                         nsim = 8,
                                         cl = cl,
                                         seed = 12))
  suppressWarnings(boot_2 <- compute_tests(fit,
                                           test_intercepts = T,
                                           test_ran_ef = T,
                                           type = 3,
                                           tests = "bootstrap",
                                           nsim = 8,
                                           #cl = cl,
                                           seed = 12))

  expect_equal(unlist(boot_1$pb_test_results[1:4, 1]), unlist(boot_2$pb_test_results[1:4, 1]))
  expect_equal(unlist(boot_1$pb_test_results[1:4, 3]), unlist(boot_2$pb_test_results[1:4, 3]))
  expect_equal(boot_1$seed, boot_2$seed)

  stopCluster(cl)
})

