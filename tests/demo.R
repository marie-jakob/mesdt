##------------------------------------------------------------
## Demo mesdt
## Date: September 2024
## Marie Jakob <marie.jakob@psychologie.uni-freiburg.de>
## ------------------------------------------------------------

library(tidyverse)
library(parallel)
library(lme4)
library(glmmTMB)

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
#### Intended Workflow ####

# 1. Fit the model
# set backend to lme4
options("mesdt.backend" = "lme4")

fit_sdt <- fit_mesdt(formula_mu = ~ emp_gender + (1 | id),
                     formula_lambda = ~ committee + (committee | id),
                     data = dat_exp_1,
                     dv = "response",
                     trial_type_var = "status_fac")
summary(fit_sdt$fit_obj)
# Model object now contains every information from the fitted model (including Wald tests)
# Summary method prints coefficients and Wald tests for the parameters


# backend glmmTMB
options("mesdt.backend" = "glmmTMB")
fit_sdt_glmmTMB <- fit_mesdt(formula_mu = ~ emp_gender * participant_gender + (1 | id),
                             formula_lambda = ~ committee * emp_gender * participant_gender + (committee | id),
                             data = dat_exp_1,
                             dv = "response",
                             trial_type_var = "status_fac")
summary(fit_sdt_glmmTMB)


# 2. Compute relevant tests
# e.g., LRTs
CL <- makeCluster(4, type = "SOCK")

fit_sdt_glmmTMB <- compute_tests(fit_sdt_glmmTMB,
                        data = dat_exp_1,
                        type = 3,
                        cl = CL,
                        test_intercepts = T,
                        test_params_mu = "all",
                        test_params_lambda = ~ committee,
                        test_ran_ef = F)
# -> Method returns the modified fit object with the added tests
# Summary method prints coefficients and all tests that are stored in the model

summary(fit_sdt_glmmTMB)

# or PB tests
fit_sdt_glmmTMB <- compute_tests(fit_sdt_glmmTMB,
                         data = dat_exp_1,
                         type = 3,
                         cl = CL,
                         tests = "bootstrap",
                         nsim = 2,
                         test_intercepts = F,
                         test_params_mu = ~ participant_gender,
                         test_params_lambda = ~ committee,
                         test_ran_ef = F)
stopCluster(CL)
summary(fit_sdt_glmmTMB)
# Summary method now also prints PB tests

# (2b. Tests for random effects -> only Type 3)

# 3. Post-processing with emmeans()

# -> Formula from


