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

  # random effects estimates
  expect_equal(ranef(fit)$ID[, 1], ranef(model_test)$ID[, 1], tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 2], ranef(model_test)$ID[, 2], tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 3], ranef(model_test)$ID[, 3] * 2, tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 4], ranef(model_test)$ID[, 4] * 2, tolerance = 1e-3)
}
)


test_that("fit_mlsdt() works for uncorrelated mu and lambda effects", {
  #fit <- fit_mlsdt(~ 1 + x1 + (1 + x1 || ID), ~ 1 + x1 + (1 + x1 || ID), dv = "y", data = internal_sdt_data)$fit_obj

  # Number of estimated fixed effects parameters
  #expect_equal(length(fixef(fit)), length(fixef(model_test_uncor)))
  # Number of estimated random effects parameters
  #expect_equal(length(ranef(fit)), length(ranef(model_test_uncor)))
  #expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(model_test_uncor))))

  # fixed effects estimates
  #expect_equal(unname(fixef(fit))[1:2], unname(fixef(model_test_uncor))[1:2], tolerance = 1e-4)
  # mu fixef effects
  #expect_equal(unname(fixef(fit))[3:4], unname(fixef(model_test_uncor))[3:4] * 2, tolerance = 1e-4)

  # observed Fisher information
  #expect_equal(unname(vcov(fit))[1:2, 1:2], unname(vcov(model_test_uncor))[1:2, 1:2], tolerance = 1e-3)
  #expect_equal(unname(vcov(fit))[3:4, 3:4], unname(vcov(model_test_uncor))[3:4, 3:4] * 4, tolerance = 1e-3)
  #expect_equal(unname(vcov(fit))[1:2, 3:4], unname(vcov(model_test_uncor))[1:2, 3:4] * 2, tolerance = 1e-3)

  # lambda random effect variances
  #expect_equal(as.data.frame(VarCorr(fit))$vcov[1:2], as.data.frame(VarCorr(model_test_uncor))$vcov[1:2], tolerance = 1e-3)
  #expect_equal(as.data.frame(VarCorr(fit))$vcov[3:4], as.data.frame(VarCorr(model_test_uncor))$vcov[3:4] * 4, tolerance = 1e-3)

  # random effects correlations
  #expect_equal(as.data.frame(VarCorr(fit))$sdcor[5:10],
  #             as.data.frame(VarCorr(model_test_uncor))$sdcor[5:10],
  #             tolerance = 1e-3)

  # random effects estimates
  #expect_equal(ranef(fit)$ID[, 1], ranef(model_test_uncor)$ID[, 1], tolerance = 1e-4)
  #expect_equal(ranef(fit)$ID[, 2], ranef(model_test_uncor)$ID[, 2], tolerance = 1e-4)
  #expect_equal(ranef(fit)$ID[, 3], ranef(model_test_uncor)$ID[, 3] * 2, tolerance = 1e-4)
  #expect_equal(ranef(fit)$ID[, 4], ranef(model_test_uncor)$ID[, 4] * 2, tolerance = 1e-4)

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

  # random effects estimates
  expect_equal(ranef(fit)$ID[, 1], ranef(model_test_uncor)$ID[, 1], tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 2], ranef(model_test_uncor)$ID[, 2], tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 3], ranef(model_test_uncor)$ID[, 3] * 2, tolerance = 1e-4)
  expect_equal(ranef(fit)$ID[, 4], ranef(model_test_uncor)$ID[, 4] * 2, tolerance = 1e-4)

}
)

test_that("fit_mlsdt() notifies the user that only uncorrelated or correlated
          random effects are possible at the moment, iff specified.", {
  expect_message(fit_mlsdt(~ 1 + x1 + (1 + x1 || ID), ~ 1 + x1 + (1 + x1 | ID), dv = "y", data = internal_sdt_data))
  expect_message(fit_mlsdt(~ 1 + x1 + (1 | ID) + (x1 | ID), ~ 1 + x1 + (1 + x1 || ID), dv = "y", data = internal_sdt_data))

  expect_no_message(fit_mlsdt(~ 1 + x1 + (x1 | ID), ~ 1 + x1 + (x1 | ID), dv = "y", data = internal_sdt_data))
  expect_no_message(fit_mlsdt(~ 1 + x1 + (1 | ID), ~ 1 + x1 + (1 | ID), dv = "y", data = internal_sdt_data))

})

#------------------------------------------------------------------------------#
#### Fit Model with Crossed Random Effects ####

test_that("fit_mlsdt() works for crossed random effects with random intercepts and no predictors", {
  fit <- fit_mlsdt(~ 1 + (1 | id) + (1 | file_name),
                   ~ 1 + (1 | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_ef")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(fit_cross_intercept)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(fit_cross_intercept)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(fit_cross_intercept))))
  # dfs & LL
  expect_equal(df.residual(fit), df.residual(fit_cross_intercept))
  expect_equal(logLik(fit), logLik(fit_cross_intercept), tolerance = 1e-3)

  # fixed effects estimates
  expect_equal(unname(fixef(fit))[1], unname(fixef(fit_cross_intercept))[1], tolerance = 1e-4)
  # mu fixef effects
  expect_equal(unname(fixef(fit))[2], unname(fixef(fit_cross_intercept))[2] * 2, tolerance = 1e-4)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1, 1], unname(vcov(fit_cross_intercept))[1, 1], tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[1, 2], unname(vcov(fit_cross_intercept))[1, 2] * 2, tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[2, 1], unname(vcov(fit_cross_intercept))[2, 1] * 2, tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[2, 2], unname(vcov(fit_cross_intercept))[2, 2] * 4, tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(as.data.frame(VarCorr(fit))$vcov[1], as.data.frame(VarCorr(fit_cross_intercept))$vcov[1] * 4, tolerance = 1e-2)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[2], as.data.frame(VarCorr(fit_cross_intercept))$vcov[2], tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[3], as.data.frame(VarCorr(fit_cross_intercept))$vcov[3] * 4, tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[4], as.data.frame(VarCorr(fit_cross_intercept))$vcov[4] * 2, tolerance = 1e-3)

  # random effects estimates
  expect_equal(ranef(fit)$file_name[, 1], ranef(fit_cross_intercept)$file_name[, 1] * 2, tolerance = 1e-2)
  expect_equal(ranef(fit)$id[, 1], ranef(fit_cross_intercept)$id[, 1], tolerance = 1e-4)
  expect_equal(ranef(fit)$id[, 2], ranef(fit_cross_intercept)$id[, 2] * 2, tolerance = 1e-4)
}
)


test_that("fit_mlsdt() works for crossed random effects with random intercepts, predictors and random slopes", {
  fit <- fit_mlsdt(~ committee + (1 | id) + (1 | file_name),
                   ~ committee + (committee | id),
                   dv = "assessment", data = dat_exp_2, trial_type_var = "status_ef")$fit_obj

  # Number of estimated fixed effects parameters
  expect_equal(length(fixef(fit)), length(fixef(fit_cross_slopes)))
  # Number of estimated random effects parameters
  expect_equal(length(ranef(fit)), length(ranef(fit_cross_slopes)))
  expect_equal(length(unlist(VarCorr(fit))), length(unlist(VarCorr(fit_cross_slopes))))
  # dfs
  expect_equal(df.residual(fit), df.residual(fit_cross_slopes))
  # log likelihoods
  expect_equal(logLik(fit), logLik(fit_cross_slopes), tolerance = 1e-4)

  # fixed effects estimates
  expect_equal(unname(fixef(fit))[1], unname(fixef(fit_cross_slopes))[1], tolerance = 1e-4)
  expect_equal(unname(fixef(fit))[2], unname(fixef(fit_cross_slopes))[2] * -1, tolerance = 1e-3)
  # mu fixef effects
  expect_equal(unname(fixef(fit))[3], unname(fixef(fit_cross_slopes))[3] * 2, tolerance = 1e-4)
  expect_equal(unname(fixef(fit))[4], unname(fixef(fit_cross_slopes))[4] * -2, tolerance = 1e-3)

  # observed Fisher information
  expect_equal(unname(vcov(fit))[1, 1], unname(vcov(fit_cross_slopes))[1, 1], tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[2, 2], unname(vcov(fit_cross_slopes))[2, 2], tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[3, 3], unname(vcov(fit_cross_slopes))[3, 3] * 4, tolerance = 1e-3)
  expect_equal(unname(vcov(fit))[4, 4], unname(vcov(fit_cross_slopes))[4, 4] * 4, tolerance = 1e-3)
  # test some other elements of the matrix as well
  expect_equal(unname(vcov(fit))[1, 2], unname(vcov(fit_cross_slopes))[1, 2] * -1, tolerance = 1e-2)
  expect_equal(unname(vcov(fit))[1, 3], unname(vcov(fit_cross_slopes))[1, 3] * 2, tolerance = 1e-2)
  expect_equal(unname(vcov(fit))[1, 4], unname(vcov(fit_cross_slopes))[1, 4] * -2, tolerance = 1e-3)

  # lambda random effect variances
  expect_equal(as.data.frame(VarCorr(fit))$vcov[1], as.data.frame(VarCorr(fit_cross_slopes))$vcov[1] * 4, tolerance = 1e-4)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[2], as.data.frame(VarCorr(fit_cross_slopes))$vcov[2], tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[3], as.data.frame(VarCorr(fit_cross_slopes))$vcov[3], tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$vcov[4], as.data.frame(VarCorr(fit_cross_slopes))$vcov[4] * 4, tolerance = 1e-3)
  expect_equal(as.data.frame(VarCorr(fit))$sdcor[6], as.data.frame(VarCorr(fit_cross_slopes))$sdcor[6], tolerance = 1e-2)
  expect_equal(as.data.frame(VarCorr(fit))$sdcor[7], as.data.frame(VarCorr(fit_cross_slopes))$sdcor[7] * -1, tolerance = 1e-1)

  # random effects estimates
  # -> does not work here for singular fits...
}
)
