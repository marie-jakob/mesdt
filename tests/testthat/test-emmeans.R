


test_that("emmeans.mesdt_fit() gives the same results as emmeans for lme4", {

  options("mesdt.backend" = "lme4")

  library(lme4)
  library(emmeans)

  # 1. One factor
  test_mod_lme <- glmer(assessment ~ status_fac * committee  +
                         (status_fac | id),
                       data = dat_exp_2,
                       family = binomial("probit"),
                       nAGQ = 0)
  test_mod_mesdt <- fit_mesdt(~ committee + (1 | id),
                              ~ committee + (1 | id),
                              data = dat_exp_2,
                              trial_type_var = "status_fac",
                              dv = "assessment",
                              correlate_sdt_params = T)
  logLik(test_mod_lme)
  logLik(test_mod_mesdt$fit_obj)
  df.residual(test_mod_lme)
  df.residual(test_mod_mesdt$fit_obj)
  fixef(test_mod_lme)
  fixef(test_mod_mesdt$fit_obj)

  # Response Bias
  emm_lme_c <- data.frame(emmeans(test_mod_lme, ~ committee))
  emm_mesdt_c <- data.frame(emmeans(test_mod_mesdt, ~ committee, dpar = "response bias"))
  expect_equal(emm_lme_c$emmean, emm_mesdt_c$emmean, tolerance = 1e-4)
  expect_equal(emm_lme_c$se, emm_mesdt_c$se, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.LCL, emm_mesdt_c$asymp.LCL, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.UCL, emm_mesdt_c$asymp.UCL, tolerance = 1e-3)

  # Discriminability
  emm_lme <- data.frame(contrast(emmeans(test_mod_lme, ~ status_fac * committee),
                      list("denied" = c(-1, 1, 0, 0),
                           "granted" = c(0, 0, -1, 1))))
  emm_mesdt <- data.frame(emmeans(test_mod_mesdt, ~ committee, dpar = "sensitivity"))
  expect_equal(emm_lme$estimate, emm_mesdt$emmean, tolerance = 1e-4)

  # 2. One factor with three levels
  test_mod_lme <- glmer(assessment ~ status_fac * contingencies +
                          (status_fac | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)

  test_mod_mesdt <- fit_mesdt(~ contingencies + (1 | id),
                              ~ contingencies + (1 | id),
                              data = dat_exp_2,
                              trial_type_var = "status_fac",
                              dv = "assessment",
                              correlate_sdt_params = T)
  # Response Bias
  emm_lme_c <- data.frame(emmeans(test_mod_lme, ~ contingencies))
  emm_mesdt_c <- data.frame(emmeans(test_mod_mesdt, ~ contingencies, dpar = "response bias"))
  expect_equal(emm_lme_c$emmean, emm_mesdt_c$emmean, tolerance = 1e-4)
  expect_equal(emm_lme_c$se, emm_mesdt_c$se, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.LCL, emm_mesdt_c$asymp.LCL, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.UCL, emm_mesdt_c$asymp.UCL, tolerance = 1e-4)

  # Sensitivity:
  emm_lme <- data.frame(contrast(emmeans(test_mod_lme, ~ status_fac * contingencies),
                      list("regular" = c(-1, 1, 0, 0, 0, 0),
                           "balanced" = c(0, 0, -1, 1, 0, 0),
                           "reversed" = c(0, 0, 0, 0, -1, 1))))
  emm_mesdt <- data.frame(emmeans(test_mod_mesdt, ~ contingencies, dpar = "sensitivity"))
  expect_equal(emm_lme$estimate, emm_mesdt$emmean, tolerance = 1e-4)
  expect_equal(emm_lme$SE, emm_mesdt$SE, tolerance = 1e-4)

  # 3. 2-Faktoriell
  test_mod_lme <- glmer(assessment ~ status_fac * contingencies * emp_gender +
                          (status_fac | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)

  test_mod_mesdt <- fit_mesdt(~ contingencies * emp_gender + (1 | id),
                              ~ contingencies * emp_gender + (1 | id),
                              data = dat_exp_2,
                              trial_type_var = "status_fac",
                              dv = "assessment",
                              correlate_sdt_params = T)
  # Response Bias
  emm_lme_c <- data.frame(emmeans(test_mod_lme, ~ contingencies * emp_gender))
  emm_mesdt_c <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "response bias"))
  expect_equal(emm_lme_c$emmean, emm_mesdt_c$emmean, tolerance = 1e-4)
  expect_equal(emm_lme_c$se, emm_mesdt_c$se, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.LCL, emm_mesdt_c$asymp.LCL, tolerance = 1e-4)
  expect_equal(emm_lme_c$asymp.UCL, emm_mesdt_c$asymp.UCL, tolerance = 1e-4)

  # Sensitivity:
  emm_lme <- data.frame(contrast(emmeans(test_mod_lme, ~ status_fac * contingencies * emp_gender),
                                 list("1" = c(-1, 1, rep(0, 10)),
                                      "2" = c(rep(0, 2), -1, 1, rep(0, 8)),
                                      "3" = c(rep(0, 4), -1, 1, rep(0, 6)),
                                      "4" = c(rep(0, 6), -1, 1, rep(0, 4)),
                                      "5" = c(rep(0, 8), -1, 1, rep(0, 2)),
                                      "6" = c(rep(0, 10), -1, 1))))
  emm_mesdt <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "sensitivity"))
  expect_equal(emm_lme$estimate, emm_mesdt$emmean, tolerance = 1e-4)
  expect_equal(emm_lme$SE, emm_mesdt$SE, tolerance = 1e-4)

})


test_that("emmeans.mesdt_fit() works for glmmTMB", {
  library(emmeans)
  library(glmmTMB)
  options("mesdt.backend" = "glmmTMB")
  test_mod_tmb <- glmmTMB(assessment ~ status_fac * contingencies * emp_gender +
                          (status_fac | id),
                        data = dat_exp_2,
                        family = binomial("probit"))

  test_mod_mesdt <- fit_mesdt(~ contingencies * emp_gender + (1 | id),
                              ~ contingencies * emp_gender + (1 | id),
                              data = dat_exp_2,
                              trial_type_var = "status_fac",
                              dv = "assessment",
                              correlate_sdt_params = T)

  # Response Bias
  emm_tmb_c <- data.frame(emmeans(test_mod_tmb, ~ contingencies * emp_gender))
  emm_mesdt_c <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "response bias"))
  expect_equal(emm_tmb_c$emmean, emm_mesdt_c$emmean, tolerance = 1e-2)
  expect_equal(emm_tmb_c$se, emm_mesdt_c$se, tolerance = 1e-4)
  expect_equal(emm_tmb_c$asymp.LCL, emm_mesdt_c$asymp.LCL, tolerance = 1e-2)
  expect_equal(emm_tmb_c$asymp.UCL, emm_mesdt_c$asymp.UCL, tolerance = 1e-2)

  # Sensitivity:
  emm_tmb <- data.frame(contrast(emmeans(test_mod_tmb, ~ status_fac * contingencies * emp_gender),
                                 list("1" = c(-1, 1, rep(0, 10)),
                                      "2" = c(rep(0, 2), -1, 1, rep(0, 8)),
                                      "3" = c(rep(0, 4), -1, 1, rep(0, 6)),
                                      "4" = c(rep(0, 6), -1, 1, rep(0, 4)),
                                      "5" = c(rep(0, 8), -1, 1, rep(0, 2)),
                                      "6" = c(rep(0, 10), -1, 1))))
  emm_mesdt <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "sensitivity"))
  expect_equal(emm_tmb$estimate, emm_mesdt$emmean, tolerance = 1e-2)
  expect_equal(emm_tmb$SE, emm_mesdt$SE, tolerance = 1e-2)
})


test_that("emmeans.mesdt_fit() works for a single-level model fit with glm()", {
  library(emmeans)
  test_mod_glm <- glmmTMB(assessment ~ status_fac * contingencies * emp_gender,
                          data = dat_exp_2,
                          family = binomial("probit"))

  test_mod_mesdt <- fit_mesdt(~ contingencies * emp_gender,
                              ~ contingencies * emp_gender,
                              data = dat_exp_2,
                              trial_type_var = "status_fac",
                              dv = "assessment",
                              correlate_sdt_params = T)

  # Response Bias
  emm_glm_c <- data.frame(emmeans(test_mod_glm, ~ contingencies * emp_gender))
  emm_mesdt_c <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "response bias"))
  expect_equal(emm_glm_c$emmean, emm_mesdt_c$emmean, tolerance = 1e-2)
  expect_equal(emm_glm_c$se, emm_mesdt_c$se, tolerance = 1e-4)
  expect_equal(emm_glm_c$asymp.LCL, emm_mesdt_c$asymp.LCL, tolerance = 1e-2)
  expect_equal(emm_glm_c$asymp.UCL, emm_mesdt_c$asymp.UCL, tolerance = 1e-2)

  # Sensitivity:
  emm_glm <- data.frame(contrast(emmeans(test_mod_glm, ~ status_fac * contingencies * emp_gender),
                                 list("1" = c(-1, 1, rep(0, 10)),
                                      "2" = c(rep(0, 2), -1, 1, rep(0, 8)),
                                      "3" = c(rep(0, 4), -1, 1, rep(0, 6)),
                                      "4" = c(rep(0, 6), -1, 1, rep(0, 4)),
                                      "5" = c(rep(0, 8), -1, 1, rep(0, 2)),
                                      "6" = c(rep(0, 10), -1, 1))))
  emm_mesdt <- data.frame(emmeans(test_mod_mesdt, ~ contingencies * emp_gender, dpar = "sensitivity"))
  expect_equal(emm_glm$estimate, emm_mesdt$emmean, tolerance = 1e-2)
  expect_equal(emm_glm$SE, emm_mesdt$SE, tolerance = 1e-2)


})

