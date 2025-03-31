##------------------------------------------------------------
## Demo mlsdt
## Date: September 2024
## Marie Jakob <marie.jakob@psychologie.uni-freiburg.de>
## ------------------------------------------------------------

library(tidyverse)
library(mlsdt)

dat_exp_1 <- readRDS("tests/test-dat/data_prep.rds") %>%
  dplyr::mutate(committee = factor(committee_ef),
                emp_gender = factor(emp_gender_ef),
                status_fac = factor(status_ef),
                status_ef = status_ef / 2,
                participant_gender = factor(vp_gender_ef))

contrasts(dat_exp_1$committee) <- contr.sum(2)
contrasts(dat_exp_1$emp_gender) <- contr.sum(2)
contrasts(dat_exp_1$participant_gender) <- contr.sum(2)
contrasts(dat_exp_1$status_fac) <- contr.sum(2)


#------------------------------------------------------------------------------#
#### Fit SDT models ####

# set backend to lme4
options("mlsdt.backend" = "lme4")

fit_sdt <- fit_mlsdt(formula_mu = ~ emp_gender * participant_gender + (1 | id),
                     formula_lambda = ~ committee * emp_gender * participant_gender + (committee | id),
                     data = dat_exp_1,
                     dv = "response",
                     trial_type_var = "status_fac")
summary(fit_sdt$fit_obj)
fit_sdt$Lambda
fit_sdt$Mu

options("mlsdt.backend" = "glmmTMB")
fit_sdt_glmmTMB <- fit_mlsdt(formula_mu = ~ emp_gender * participant_gender + (1 | id),
                     formula_lambda = ~ committee * emp_gender * participant_gender + (committee | id),
                     data = dat_exp_1,
                     dv = "response",
                     trial_type_var = "status_fac")


LRTs <- compute_tests(fit_sdt_glmmTMB,
                      data = dat_exp_1,
                      type = 3,
                      test_intercepts = T,
                      test_params_mu = "all",
                      test_params_lambda = ~ committee,
                      test_ran_ef = F)

LRTs$LRTs


# Test random effects
LRTs_ranef <- compute_tests(fit_sdt_glmmTMB,
                            data = dat_exp_1,
                            type = 3,
                            test_intercepts = T,
                            test_ran_ef = F)

# Significance tests based on parametric bootstrapping
boot_tests <- compute_tests(fit_sdt_glmmTMB,
                      data = dat_exp_1,
                      tests = "bootstrap",
                      nsim = 100,
                      type = 3,
                      test_intercepts = T,
                      test_ran_ef = F)


