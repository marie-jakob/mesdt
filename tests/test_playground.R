

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



dat_exp_2 <- readRDS("tests/test-dat/dat_exp_2.rds") %>%
  dplyr::filter(id %in% c("regular-1", "regular-2",
                   "reversed-1", "reversed-2",
                   "balanced-1", "balanced-2"))

# Use a reduced model that only includes two terms

fit_exp_2 <- fit_mlsdt(~ committee_ef * emp_gender_ef + (1 | id),
                       ~ committee_ef * emp_gender_ef + (1 | id),
                       dv = "assessment",
                       trial_type_var = "status_ef",
                       data = dat_exp_2)

mm_exp_2 <- construct_modelmatrices(~ committee_ef * emp_gender_ef + (1 | id),
                                    ~ committee_ef * emp_gender_ef + (1 | id),
                                    dv = "assessment",
                                    trial_type_var = "status_ef",
                                    data = dat_exp_2)

LRTs_exp_2_2 <- compute_LRTs(fit_exp_2$fit_obj,
                           ~ committee_ef * emp_gender_ef + (1 | id),
                           ~ committee_ef * emp_gender_ef + (1 | id),
                           dv = "assessment",
                           data = dat_exp_2,
                           type = 2,
                           mm_exp_2)

LRTs_exp_2_3 <- compute_LRTs(fit_exp_2$fit_obj,
                           ~ committee_ef * emp_gender_ef + (1 | id),
                           ~ committee_ef * emp_gender_ef + (1 | id),
                           dv = "assessment",
                           data = dat_exp_2,
                           type = 3,
                           mm_exp_2)

summary(fit_exp_2$fit_obj)

LRTs_exp_2_2$LRTs
LRTs_exp_2_3$LRTs
