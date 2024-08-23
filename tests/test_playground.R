

test_formula <- construct_glmer_formula(
  formula_mu = ~ 1 + x1 + (x1 || ID),
  formula_lambda = ~ 1 + x1 + (x1 | ID),
  dv = "dv",
  mm = construct_modelmatrices(~ 1 + x1 + (x1 || ID), ~ 1 + x1 + (x1 || ID), data = internal_sdt_data),
)


dv <- "y"

mm <- construct_modelmatrices(~ 1 + x1 + (x1 | ID),
                              ~ 1 + x1 + (x1 | ID),
                              data = internal_sdt_data)


fit <- fit_mlsdt(~ x1 + (x1 | ID), ~ x1 + (x1 | ID), dv = "y", data = internal_sdt_data)

summary(fit$fit_obj)


new_formula <- reduce_random_effects_structure(test_formula, fit$fit_obj)



fit_2 <- lme4::glmer(test_formula,
                     data = internal_sdt_data,
                     family = binomial(link = "probit"),
                     # this is only for testing speed -> changed for actual use
                     nAGQ = 0)


fit$Mu

fit$fit_obj



dat_rev <- readRDS("tests/test-dat/data_rev.rds")


fit_test <- fit_mlsdt(~ emp_gender_ef * participant_gender_ef * condition_ef +
                        (1 | id),
                      ~ committee_ef * emp_gender_ef * participant_gender_ef * condition_ef +
                        (committee_ef || id),
                      dv = "response_dum",
                      trial_type_var = "status",
                      data = dat_rev)
fit_test$Mu
fit_test$Lambda


summary(fit_test$fit_obj)




#------------------------------------------------------------------------------#
#### Test type 2 LRTs ####


# Use a reduced model that only includes two terms

form_lambda <- ~ committee + (1 | id)
form_mu <- ~ committee * emp_gender + (1 | id)

fit_exp_2 <- fit_mlsdt(form_mu,
                       form_lambda,
                       dv = "assessment",
                       trial_type_var = "status",
                       data = dat_exp_2)

mm_exp_2 <- construct_modelmatrices(form_mu,
                                    form_lambda,
                                    dv = "assessment",
                                    trial_type_var = "status",
                                    data = dat_exp_2)

LRTs_exp_2_2 <- compute_LRTs(fit_exp_2$fit_obj,
                             form_mu,
                             form_lambda,
                             dv = "assessment",
                             data = dat_exp_2,
                             type = 2,
                             mm_exp_2)


# Compute LRTs manually for testing
# full model - main effects lambda
full_model_main_effects <- glmer(assessment ~ emp_gender_ef + status_ef + status_ef:committee_ef +
                                   status_ef:emp_gender_ef + status_ef:emp_gender_ef:committee_ef + (status_ef | id),
                                 data = dat_exp_2,
                                 family = binomial("probit"),
                                 nAGQ = 0)
# lambda - committee_ef
lambda_committee_red <- glmer(assessment ~ status_ef + status_ef:committee_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                              data = dat_exp_2,
                              family = binomial("probit"),
                              nAGQ = 0)
# lambda - emp_gender_ef
lambda_emp_gender_red <- glmer(assessment ~ committee_ef * status_ef + status_ef:emp_gender_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                               data = dat_exp_2,
                               family = binomial("probit"),
                               nAGQ = 0)

# lambda - emp_gender_ef:committee_ef
# full model
model_full <- glmer(assessment ~ committee_ef * emp_gender_ef * status_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)
lambda_interaction <- glmer(assessment ~ committee_ef * status_ef + emp_gender_ef * status_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)

#### Mu
# full model - main effects lambda
full_model_main_effects_mu <- glmer(assessment ~ emp_gender_ef * committee_ef + status_ef + status_ef:committee_ef +
                                      status_ef:emp_gender_ef + (status_ef | id),
                                    data = dat_exp_2,
                                    family = binomial("probit"),
                                    nAGQ = 0)
# lambda - committee_ef
mu_committee_red <- glmer(assessment ~ emp_gender_ef * committee_ef + status_ef * emp_gender_ef + (status_ef | id),
                          data = dat_exp_2,
                          family = binomial("probit"),
                          nAGQ = 0)

# lambda - emp_gender_ef
mu_emp_gender_red <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef * committee_ef + (status_ef | id),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)

# lambda - emp_gender_ef:committee_ef
# full model
mu_interaction <- glmer(assessment ~ committee_ef * emp_gender_ef + emp_gender_ef * status_ef + committee_ef * status_ef + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)


anova(lambda_committee_red, full_model_main_effects)
anova(lambda_emp_gender_red, full_model_main_effects)
anova(model_full, lambda_interaction)
anova(mu_committee_red, full_model_main_effects_mu)
anova(mu_emp_gender_red, full_model_main_effects_mu)
anova(model_full, mu_interaction)

chisquares_two_factors <- c(
  -2 * (logLik(lambda_committee_red) - logLik(full_model_main_effects)),
  -2 * (logLik(lambda_emp_gender_red) - logLik(full_model_main_effects)),
  -2 * (logLik(model_full) - logLik(lambda_interaction)),
  -2 * (logLik(mu_committee_red) - logLik(full_model_main_effects_mu)),
  -2 * (logLik(mu_emp_gender_red) - logLik(full_model_main_effects_mu)),
  -2 * (logLik(model_full) - logLik(mu_interaction))
)


#------------------------------------------------------------------------------#
#### Type II LRTs only main effects ####

# -> should be identical for Type II and Type III SS
# and identical to afex Type III
form_lambda <- ~ emp_gender * committee + (1 | id)
form_mu <- ~ emp_gender * committee + (1 | id)

fit_exp_2 <- fit_mlsdt(form_mu,
                       form_lambda,
                       dv = "assessment",
                       trial_type_var = "status",
                       data = dat_exp_2)

mm_exp_2 <- construct_modelmatrices(form_mu,
                                    form_lambda,
                                    dv = "assessment",
                                    trial_type_var = "status",
                                    data = dat_exp_2)



LRTs_exp_2_2 <- compute_LRTs(fit_exp_2$fit_obj,
                             form_mu,
                             form_lambda,
                             dv = "assessment",
                             data = dat_exp_2,
                             type = 2,
                             mm_exp_2,
                             test_intercepts = F)


LRTs_exp_2_3 <- compute_LRTs(fit_exp_2$fit_obj,
                             form_mu,
                             form_lambda,
                             dv = "assessment",
                             data = dat_exp_2,
                             type = 3,
                             mm_exp_2,
                             test_intercepts = T)

LRTs_exp_2_3
fit_exp_2_afex
LRTs_exp_2_2


# Tests:
# Only intercepts
# one intercept + one predictor
# one predictor mu + lambda


#------------------------------------------------------------------------------#
#### Factors with > 2 levels ####

form_lambda <- ~ 1 + (1 | id)
form_mu <- ~ 1 + (1 | id)

fit_exp_2 <- fit_mlsdt(form_mu,
                       form_lambda,
                       dv = "assessment",
                       trial_type_var = "status",
                       data = dat_exp_2)

test <- glmer(assessment ~ status_ef + (status_ef | id),
              data = dat_exp_2,
              family = binomial("probit"),
              nAGQ = 0)

ranef(fit_exp_2$fit_obj)
ranef(test)

mm_exp_2 <- construct_modelmatrices(form_mu,
                                    form_lambda,
                                    dv = "assessment",
                                    trial_type_var = "status",
                                    data = dat_exp_2)

LRTs_exp_2_2 <- compute_LRTs(fit_exp_2$fit_obj,
                             form_mu,
                             form_lambda,
                             dv = "assessment",
                             data = dat_exp_2,
                             type = 3,
                             mm_exp_2,
                             test_intercepts = T)

# This somehow works??


# Goal for LRTs: option to test factors and to test single predictors
# -> Sometime in the future

#------------------------------------------------------------------------------#
#### Crossed random effects ####

form_cross <- assessment ~ status_ef + (1 | id) + (status_ef | id) + (1 | stimulus)
lme4::findbars(form_cross)
rdm_pred_lambda <- paste(sapply(lme4::findbars(form_cross), function(x) {
  gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
}), collapse = "+")

# find a way to split random effects according to random-effects grouping factor

mod_cross <- glmer(assessment ~ status_ef + (1 | id) + (1 | stimulus),
                   data = dat_exp_2,
                   family = binomial("probit"),
                   nAGQ = 0)

summary(mod_cross)

rdm_pred_lambda <- paste(sapply(lme4::findbars(form_cross), function(x) {
  gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
}), collapse = "+")

mm_rdm_lambda <- stats::model.matrix(formula(paste("~", rdm_pred_lambda, sep = "")),
                                     data = dat_exp_2)


