##------------------------------------------------------------
## Demo mlsdt
## Date: September 2024
## Marie Jakob <marie.jakob@psychologie.uni-freiburg.de>
## ------------------------------------------------------------

library(tidyverse)

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


LRTs <- compute_LRTs(fit_obj = fit_sdt_glmmTMB$fit_obj,
                     mm = construct_modelmatrices(
                       formula_mu = ~ emp_gender * participant_gender + (1 | id),
                       formula_lambda = ~ committee * emp_gender * participant_gender + (committee | id),
                       data = dat_exp_1,
                       trial_type_var = "status_fac",
                     ),
                     formula_mu = ~ emp_gender * participant_gender + (1 | id),
                     formula_lambda = ~ committee * emp_gender * participant_gender + (committee | id),
                     data = dat_exp_1,
                     # choose between type II and type III sums of squares
                     type = 3,
                     dv = "response",
                     # optionally test intercepts for sensitivity and response bias
                     test_intercepts = T)

LRTs$LRTs

# Test random effects (done, but not yet committed)

# Automatically reduced random-effects structure (two strategies: maximal & parsimonious)

# Significance tests based on parametric bootstrapping

# User-friendly helper functions


