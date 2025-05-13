
options("mesdt.backend" = "lme4")

#------------------------------------------------------------------------------#
#### fit_mesdt() ####

# use a saved model for this

test_that("fit_mesdt() estimates the correct model", {
  fit <- fit_mesdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data,
                   trial_type_var = "trial_type_fac")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(model_test)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(model_test)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test))))

  # dfs & LL
  expect_equal(df.residual(fit), df.residual(model_test))
  expect_equal(logLik(fit), logLik(model_test), tolerance = 1e-3)

  # fixed effects estimates
  expect_equal(unname(fixef(fit)), unname(fixef(model_test)), tolerance = 1e-6)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1:4, 1:4], unname(vcov(model_test))[1:4, 1:4], tolerance = 1e-1)

  # random effect variances and covariance
  expect_equal(as.data.frame(VarCorr(fit))$vcov, as.data.frame(VarCorr(model_test))$vcov, tolerance = 1e-5)

  # random effects estimates
  expect_equal(unname(ranef(fit)$ID), unname(ranef(model_test)$ID), tolerance = 1e-4)
}
)


test_that("fit_mesdt() works for uncorrelated mu and lambda effects", {
  fit <- fit_mesdt(~ 1 + committee + (1 + committee | id),
                   ~ 1 + committee + (1 + committee | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac",
                   correlate_sdt_params = F)$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(model_uncor_sdt)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(model_uncor_sdt)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_uncor_sdt))))

  # dfs & LL
  expect_equal(df.residual(fit), df.residual(model_uncor_sdt))
  expect_equal(logLik(fit), logLik(model_uncor_sdt), tolerance = 1e-3)

  # fixed effects estimates
  expect_equal(abs(unname(fixef(fit))), abs(unname(fixef(model_uncor_sdt))), tolerance = 1e-4)

  # observed Fisher information
  expect_equal(abs(unname(vcov(fit))[1:4, 1:4]), abs(unname(vcov(model_uncor_sdt))[1:4, 1:4]), tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(abs(as.data.frame(VarCorr(fit))$vcov), abs(as.data.frame(VarCorr(model_uncor_sdt))$vcov), tolerance = 1e-3)

  # random effects estimates
  expect_equal(unname(abs(ranef(fit)$id[, 1:4])), unname(abs(ranef(model_uncor_sdt)$id[, 1:4])), tolerance = 1e-4)

}
)


test_that("fit_mesdt() works for uncorrelated random effects (|| notation)", {
  fit <- fit_mesdt(~ 1 + x1 + (1 + x1 || ID),
                   ~ 1 + x1 + (1 + x1 || ID), dv = "y",
                   data = internal_sdt_data,
                   trial_type_var = "trial_type_fac")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(model_test_uncor)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(model_test_uncor)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test_uncor))))

  # dfs & LL
  expect_equal(df.residual(fit), df.residual(model_test_uncor))
  expect_equal(logLik(fit), logLik(model_test_uncor), tolerance = 1e-3)

  # fixed effects estimates
  expect_equal(unname(fixef(fit)), unname(fixef(model_test_uncor)), tolerance = 1e-6)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1:4, 1:4], unname(vcov(model_test_uncor))[1:4, 1:4], tolerance = 1e-1)

  # random effect variances and covariance
  expect_equal(as.data.frame(VarCorr(fit))$vcov, as.data.frame(VarCorr(model_test_uncor))$vcov, tolerance = 1e-5)

  # random effects estimates
  expect_equal(unname(abs(ranef(fit)$ID)[, 1:4]), unname(abs(ranef(model_test_uncor)$ID)[, 1:4]), tolerance = 1e-3)

}
)

test_that("fit_mesdt() notifies the user that only uncorrelated or correlated
          random effects are possible at the moment, iff specified.", {
  expect_message(fit_mesdt(~ 1 + x1 + (1 + x1 || ID), ~ 1 + x1 + (1 + x1 | ID), dv = "y",
                           trial_type_var = "trial_type_fac", data = internal_sdt_data))
  expect_message(fit_mesdt(~ 1 + x1 + (1 | ID) + (x1 | ID), ~ 1 + x1 + (1 + x1 || ID), dv = "y",
                           trial_type_var = "trial_type_fac", data = internal_sdt_data))

})

#------------------------------------------------------------------------------#
#### Fit Model with Crossed Random Effects ####

test_that("fit_mesdt() works for crossed random effects with random intercepts and no predictors", {
  fit <- fit_mesdt(~ 1 + (1 | id) + (1 | file_name),
                   ~ 1 + (1 | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(fit_cross_intercept)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(fit_cross_intercept)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(fit_cross_intercept))))
  # dfs & LL
  expect_equal(df.residual(fit), df.residual(fit_cross_intercept))
  expect_equal(logLik(fit), logLik(fit_cross_intercept), tolerance = 1e-3)

  # fixed effects estimates
  expect_equal(abs(unname(fixef(fit))), abs(unname(fixef(fit_cross_intercept))), tolerance = 1e-4)
  # mu fixef effects

  # observed Fisher information
  expect_equal(abs(unname(vcov(fit))[1:2, 1:2]), abs(unname(vcov(fit_cross_intercept))[1:2, 1:2]), tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(abs(as.data.frame(VarCorr(fit))$vcov)[1:4], abs(as.data.frame(VarCorr(fit_cross_intercept))$vcov[1:4]), tolerance = 1e-2)

  # random effects estimates
  expect_equal(ranef(fit)$file_name[, 1], ranef(fit_cross_intercept)$file_name[, 1], tolerance = 1e-2)
  expect_equal(unname(abs(ranef(fit)$id)[, 1:2]), unname(abs(ranef(fit_cross_intercept)$id)[, 1:2]), tolerance = 1e-4)
}
)


test_that("fit_mesdt() works for crossed random effects with random intercepts, predictors and random slopes", {
  options("mesdt.backend" = "glmmTMB")
  fit <- fit_mesdt(~ committee + (1 | id) + (1 | file_name),
                   ~ committee + (committee | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)[[1]]), length(fixef(fit_cross_slopes)[[1]]))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)[[1]]), length(ranef(fit_cross_slopes)[[1]]))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(fit_cross_slopes))))
  # dfs
  expect_equal(df.residual(fit), df.residual(fit_cross_slopes))
  # log likelihoods
  expect_equal(logLik(fit), logLik(fit_cross_slopes), tolerance = 1e-6)

  # fixed effects estimates
  expect_equal(abs(unname(fixef(fit)[[1]])), abs(unname(fixef(fit_cross_slopes)[[1]])), tolerance = 1e-6)

  # observed Fisher information
  expect_equal(abs(unname(vcov(fit)[[1]])[1:4, 1:4]), abs(unname(vcov(fit_cross_slopes)[[1]])[1:4, 1:4]), tolerance = 1e-4)

  # lambda random effect variances
  expect_equal(abs(unlist(VarCorr(fit))), abs(unlist(VarCorr(fit_cross_slopes))), tolerance = 1e-4)
}
)


test_that("fit_mesdt() works when only mu or lambda have random effects", {
  options("mesdt.backend" = "lme4")
  fit <- fit_mesdt(~ committee,
                   ~ committee + (1 | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac")$fit_obj
  fit_test <- lme4::glmer(assessment ~ status_ef * committee + (1 | id), data = dat_exp_2, family = binomial("probit"),
                    nAGQ = 0)
  expect_equal(logLik(fit), logLik(fit_test))
  expect_equal(sort(abs(as.numeric(fixef(fit)))), sort(abs(as.numeric(fixef(fit_test)))))
  expect_equal(as.numeric(ranef(fit)), as.numeric(ranef(fit_test)))
  expect_equal(as.numeric(VarCorr(fit)), as.numeric(VarCorr(fit_test)), tolerance = 1e-5)


  fit <- fit_mesdt(~ committee + (1 | id),
                   ~ committee,
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac")$fit_obj
  fit_test <- lme4::glmer(assessment ~ status_ef * committee + (0 + status_ef | id), data = dat_exp_2, family = binomial("probit"),
                    nAGQ = 0)
  expect_equal(logLik(fit), logLik(fit_test))
  expect_equal(sort(abs(as.numeric(fixef(fit)))), sort(abs(as.numeric(fixef(fit_test)))))
  expect_equal(as.numeric(ranef(fit)), as.numeric(ranef(fit_test)))
  expect_equal(as.numeric(VarCorr(fit)), as.numeric(VarCorr(fit_test)), tolerance = 1e-5)


}
)

test_that("fit_mesdt() works without any random effects", {
  expect_message(fit_mesdt(~ committee,
                           ~ emp_gender,
                           dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac"))
  fit <- fit_mesdt(~ committee,
                   ~ committee,
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac")$fit_obj
  fit_test <- glm(assessment ~ status_ef * committee, data = dat_exp_2, family = binomial("probit"))
  expect_equal(logLik(fit), logLik(fit_test))
  expect_equal(sort(abs(as.numeric(coefficients(fit)))), sort(abs(as.numeric(coefficients(fit_test)))))
}
)


#------------------------------------------------------------------------------#
#### Test control argument ####


test_that("Setting control arguments in fit_mesdt() works", {
  options("mesdt.backend" = "lme4")
  expect_message(fit_mesdt(~ committee + (1 | id),
                           ~ emp_gender,
                           dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac",
                           control = glmmTMB::glmmTMBControl(list(iter.max=5000, eval.max=5000))))
  fit <- fit_mesdt(~ committee + (1 | id),
                   ~ emp_gender,
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac",
                   control = lme4::glmerControl(optCtrl=list(maxfun=450)))
  expect_equal(fit$fit_obj@optinfo$control$maxfun, 450)


  # Same for glmmTMB
  options("mesdt.backend" = "glmmTMB")
  expect_message(fit_mesdt(~ committee + (1 | id),
                           ~ emp_gender,
                           dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac",
                           control = lme4::glmerControl(optCtrl = list(maxfun = 500))))
  fit <- fit_mesdt(~ committee + (1 | id),
                   ~ emp_gender,
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_fac",
                   control = glmmTMB::glmmTMBControl(list(iter.max=5000)))
  expect_true("control" %in% as.character(fit$fit_obj$call))
}
)


#------------------------------------------------------------------------------#
#### Test check_sensitivity ####

test_that("mesdt throws a message when sensitivity < 0", {
  options("mesdt.backend" = "glmmTMB")
  dat_exp_2$status_rev <- ifelse(dat_exp_2$status_fac == 1, -1, 1)
  expect_warning(fit_mesdt(~ committee + (1 | id),
                         ~ emp_gender,
                         dv = "assessment", data = dat_exp_2, trial_type_var = "status_rev",
                         ))
  }
)
