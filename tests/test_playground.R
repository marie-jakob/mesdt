
library(pbkrtest)
library(lme4)
library(glmmTMB)
library(emmeans)

test_formula <- construct_glmer_formula(
  formula_mu = ~ 1 + x1 + (x1 || ID),
  formula_lambda = ~ 1 + x1 + (x1 | ID),
  dv = "dv",
  mm = construct_modelmatrices(~ 1 + x1 + (x1 || ID), ~ 1 + x1 + (x1 || ID), data = internal_sdt_data)[["mm"]],
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

form_cross <- assessment ~ status_ef + (1 | id) + (status_ef + committee_ef | id) + (1 | stimulus)
lme4::findbars(form_cross)

rdm_facs_lambda <- sapply(lme4::findbars(form_cross), function(x) {
  gsub("0 \\+", "", strsplit(as.character(x), "\\|"))[3]
})

rdm_pred_lambda <- sapply(lme4::findbars(form_cross), function(x) {
  gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
})





# find a way to split random effects according to random-effects grouping factor

mod_cross <- glmer(assessment ~ status_ef + (1 || id),
                   data = dat_exp_2,
                   family = binomial("probit"),
                   nAGQ = 0)

summary(mod_cross)

rdm_pred_lambda <- paste(sapply(lme4::findbars(form_cross), function(x) {
  gsub("0 \\+", "", strsplit(as.character(x), "\\|")[2])
}), collapse = "+")

mm_rdm_lambda <- stats::model.matrix(formula(paste("~", rdm_pred_lambda, sep = "")),
                                     data = dat_exp_2)

construct_glmer_formula(~ x1 + (1 | ID) + (1 | stim), ~ x1 + (1 | ID) + (1 | stim), dv = "trial_type",
                        correlate_sdt_params = F)


construct_modelmatrices(~ committee + (1 | id) + (1 | stimulus), ~ committee + (1 | id) + (1 | stimulus),
                        trial_type_var = "assessment", data = dat_exp_2) -> test_mm


test_cross <- fit_mlsdt(~ committee + (1 | id) + (1 | stimulus), ~ committee + (1 | id) + (1 | stimulus),
                        dv = "assessment", trial_type_var = "status", data = dat_exp_2)


#------------------------------------------------------------------------------#
#### Fit crossed random effects model ####


fit_cross_intercept <- glmer(assessment ~ status_ef + (status_ef | file_name),
                             data = dat_exp_2,
                             family = binomial("probit"),
                             nAGQ = 0)

summary(fit_cross_intercept)


test <- fit_mlsdt(~ 1 + (1 | file_name),
          ~ 1 + (1 | file_name),
          data = dat_exp_2,
          dv = "assessment",
          trial_type_var = "status")
summary(test$fit_obj)
logLik(test$fit_obj)
logLik(fit_cross_intercept)
df.residual(test$fit_obj)
df.residual(fit_cross_intercept)

ranef(fit_cross_intercept)



#------------------------------------------------------------------------------#
#### LRTs for random effects ####


# Does anova() do this correctly?
form_lambda <- ~ emp_gender * committee + (committee || id)
form_mu <- ~ emp_gender * committee + (1 | id)

fit_exp_2 <- fit_mlsdt(form_mu,
                       form_lambda,
                       dv = "assessment",
                       trial_type_var = "status_fac",
                       data = dat_exp_2)

fit_exp_2_red <- fit_mlsdt(form_lambda,
                           ~ emp_gender * committee + (1 | id),
                           dv = "assessment",
                           trial_type_var = "status_fac",
                       data = dat_exp_2)
logLik(fit_exp_2$fit_obj)
logLik(fit_exp_2_red$fit_obj)

anova(fit_exp_2$fit_obj, fit_exp_2_red$fit_obj)

# Does drop1() do this correctly?
mu_committee_red <- glmer(assessment ~ emp_gender_ef * committee_ef + status_ef * emp_gender_ef + (status_ef || id),
                          data = dat_exp_2,
                          family = binomial("probit"),
                          nAGQ = 0)
drop1(mu_committee_red, ~ (status_ef | id), test = "Chisq")
# -> No: no Chisq and p values


construct_glmer_formula(
  formula_mu = ~ 1 + x1 + (x1 | VP),
  formula_lambda = ~ 1 + x2 + (x2 | VP),
  dv = "dv",
  param_idc = 3,
  remove_from_mu = F,
  remove_from_rdm = "VP"
)


as.character(as.formula("dv ~ 0 + mm[['lambda']][, -3] + mm[['mu']] +
                            (0 + mm[['rdm_lambda_VP']] + mm[['rdm_mu_VP']] | VP)")
)



#------------------------------------------------------------------------------#
#### Parametric Bootstrapping ####




#------------------------------------------------------------------------------#
#### Random Effects Selection ####

# There are so many different opinions on this...
# -> In the manuscript: emphasize that the algorithmic implementations are _heuristic_

options("mlsdt.backend" = "lme4")

model_1 <- fit_mlsdt(formula_mu = ~ emp_gender +
                          (emp_gender | id) + (committee | file_name),
                        formula_lambda = ~ emp_gender * committee +
                          (emp_gender * committee | id) + (committee | file_name),
                        data = dat_exp_2,
                        dv = "assessment",
                        trial_type_var = "status_fac")
summary(model_1$fit_obj)
isSingular(model_1$fit_obj)
# singular fit -> remove correlations

correlations_in_model <- ifelse(any(! is.na(data.frame(VarCorr(model_1$fit_obj))$var2)), T, F)

# TODO: remove all correlations right away or (if multiple random effects grouping
# factors are present) for one factor at a time
# -> so far: only possible for all random effects at the same time
model_2 <- fit_mlsdt(formula_mu = ~ emp_gender +
                       (emp_gender || id) + (committee || file_name),
                     formula_lambda = ~ emp_gender * committee +
                       (emp_gender * committee || id) + (committee || file_name),
                     data = dat_exp_2,
                     dv = "assessment",
                     trial_type_var = "status_fac")

mm_test <- construct_modelmatrices(formula_mu = ~ emp_gender +
                                     (emp_gender || id) + (committee || file_name),
                                   formula_lambda = ~ emp_gender * committee +
                                     (emp_gender * committee || id) + (committee || file_name),
                                   data = dat_exp_2,
                                   dv = "assessment",
                                   trial_type_var = "status_fac")


variance_comps <- data.frame(VarCorr(model_2$fit_obj))
variance_comps$name <- sapply(strsplit(variance_comps$var1, '"'), function(x) return(x[2]))
variance_comps$idx <- as.numeric(sapply(variance_comps$var1, function(x) return(substr(x, nchar(x) - 2, nchar(x) -1))))



test_form <- ~ x1 * x2 + (x1 * x2 | id)


# Iterate through all random model matrices
# select the to-be-removed variance component from each
# (smallest component that can be removed while not violating marginality)
# select the smallest variance component among all selected
# generate a new formula



isSingular(model_2$fit_obj)
# singular fit
correlations_in_model <- ifelse(any(! is.na(data.frame(VarCorr(model_2$fit_obj))$var2)), T, F)

summary(model_2$fit_obj)
# Remove all variance terms that are zero
# or if there are none, the smallest ones
to_remove <- list()
var_corr <- data.frame(VarCorr(model_2$fit_obj))
idc <- which(var_corr$vcov == min(var_corr$vcov))
to_remove_idc <- sapply(strsplit(var_corr$var1[idc], "\\ "), function(x) return(as.numeric(substr(x[2], 1, 1))))

to_remove_rdm_fac <- sapply(strsplit(var_corr$grp[idc], "\\."), function(x) return(x[1]))
to_remove_sdt_param <- sapply(strsplit(var_corr$var1[idc], "\\_"), function(x) return(x[2]))

names_tmp <- paste("rdm", to_remove_sdt_param, to_remove_rdm_fac, sep = "_")

for (i in 1:length(names_tmp)) {
  if (is.null(to_remove[[names_tmp[i]]])) to_remove[[names_tmp[i]]] <- to_remove_idc[i]
  else to_remove[[names_tmp[i]]] <- c(to_remove[[names_tmp[i]]], to_remove_idc[i])
}




# Remove indices (filename intercept for mu and lambda)

model_3 <- fit_mlsdt(formula_mu = ~ emp_gender +
                       (emp_gender || id) + (0 + committee || file_name),
                     formula_lambda = ~ emp_gender * committee +
                       (emp_gender * committee || id) + (0 + committee || file_name),
                     data = dat_exp_2,
                     dv = "assessment",
                     trial_type_var = "status_fac")


summary(isSingular(model_3$fit_obj))
correlations_in_model <- ifelse(any(! is.na(data.frame(VarCorr(model_3$fit_obj))$var2)), T, F)


model_full$fit_obj$fit$message
model_full$fit_obj$fit$convergence


model_red_1$fit_obj$fit$message
model_red_1$fit_obj$fit$convergence


which(data.frame(VarCorr(model_red_1))$sdcor == min(data.frame(VarCorr(model_red_1))$sdcor))



# Algorithm maximal:

# Model: maximal model without correlations
# Repeat until convergence
# Fit model
# If model is not converged, remove variance component (according to the principle of marginality)

#------------------------------------------------------------------------------#
#### Parametric Bootstrap ####

model_full <- glmer(assessment ~ committee_ef * emp_gender_ef * status_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)

lambda_interaction <- glmer(assessment ~ committee_ef * status_ef + emp_gender_ef * status_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)

fit_cross_slopes <- glmer(assessment ~ status_ef + (status_ef | id),
                            data = dat_exp_2, family = binomial("probit"),
                          nAGQ = 0)

cross_lambda_intercept <- glmer(assessment ~ status_ef + (1 | id),
                                  data = dat_exp_2, family = binomial("probit"),
                                nAGQ = 0)

pb <- PBmodcomp(fit_cross_slopes, cross_lambda_intercept, nsim = 100)


# bootstrap with glmmTMB
model_full <- glmmTMB(assessment ~ status_ef + committee_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"))

model_red <- glmmTMB(assessment ~ status_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"))

boot_glmmTMB <- compute_parametric_bootstrap_test(model_full, model_red, nsim = 1)


# apply version

LRs_boot <- sapply(1:nsim, function(x) {
  print(x)
  sim_dat_tmp <- stats::simulate(model_red)[[1]]
  sim_fit_full <- refit(model_full, sim_dat_tmp)
  sim_fit_red <- refit(model_red, sim_dat_tmp)
  return(-2 * (logLik(sim_fit_red) - logLik(sim_fit_full)))
})


LR_emp <- -2 * (logLik(model_red) - logLik(model_full))
# according to Halekoh & Hojsgaard (2014)
p_boot <- (length(which(LRs_boot > LR_emp)) + 1) / (nsim + 1)

nsim <- 500
set.seed(3405)
start_man <- Sys.time()
boot_man <- PBmodcomp_glmmTMB(model_full, model_red, nsim = nsim)
end_man <- Sys.time()

start_pack <- Sys.time()
boot_pack <- PBmodcomp(model_full, model_red, nsim = nsim, seed = 1)
end_pack <- Sys.time()

# Time is very similar, package is probably not needed


b_full <- lme4::bootMer(model_full, FUN=function(x) -2 * logLik(x), nsim = 20, .progress="txt", type = "parametric")
b_red <- lme4::bootMer(model_red, FUN=function(x) -2 * logLik(x), nsim = 20, .progress="txt", type = "parametric")

#------------------------------------------------------------------------------#
#### PB manual ####

options("mlsdt.backend" = "lme4")
model_full <- glmmTMB(assessment ~ status_ef + committee_ef + (status_ef | id),
                      data = dat_exp_2,
                      family = binomial("probit"))

model_red <- glmmTMB(assessment ~ status_ef + (status_ef | id),
                     data = dat_exp_2,
                     family = binomial("probit"))

fit <- fit_mlsdt(~ committee_ef + (1 | id), ~ committee_ef + (1 | id),
                 data = dat_exp_2,
                 trial_type_var = "status_fac",
                 dv = "assessment")
fit_red <- fit_mlsdt(~ 1 + (1 | id), ~ committee_ef + (1 | id),
                 data = dat_exp_2,
                 trial_type_var = "status_fac",
                 dv = "assessment")


nsim <- 500
seed <- 251433
start_pack <- Sys.time()
boot_pack <- PBmodcomp(fit$fit_obj, fit_red$fit_obj, nsim = nsim, seed = seed)
end_pack <- Sys.time()



start_man <- Sys.time()
boot_man <- compute_parametric_bootstrap_test(fit$fit_obj, fit_red$fit_obj, nsim = nsim, seed = seed,
                                              mm = construct_model_matrices(~ committee_ef + (1 | id), ~ committee_ef + (1 | id),
                                                                            trial_type_var = "status_fac", dv = "assessment"),
                                              dv = "assessment",
                                              data = dat_exp_2)
end_man <- Sys.time()


#------------------------------------------------------------------------------#
#### Parallelization ####

library(parallel)
library(afex)

detectCores()
nc <- 8
#cls <- makeCluster(nc, type = "SOCK")
#clusterEvalQ(cls, {
#  library(lme4); library(afex)
#}


mixed(assessment ~ status_ef * committee_ef * emp_gender_ef + (status_ef | id),
      data = dat_exp_2,
      family = binomial("probit"),
      args_test = list(cl = cls),
      test_intercept = T,
      method = "LRT")

#------------------------------------------------------------------------------#
#### Prepare for emmeans and plots ####

dat_test <- dat_exp_2
dat_test
dat_test$committee_s <- as.numeric(as.character(dat_test$committee)) * as.numeric(as.character(dat_test$status_fac))

mod_old <- glmmTMB(assessment ~ status_fac * committee + (status_fac | id),
                      data = dat_exp_2,
                      family = binomial("probit"))

mod_new <- glmmTMB(assessment ~ status_fac + committee + committee_s + (status_fac | id),
                   data = dat_test,
                   family = binomial("probit"))
logLik(mod_old)
logLik(mod_new)
fixef(mod_old)
fixef(mod_new)
summary(mod_old)
summary(mod_new)
# same models

# emmeans for committee on sensitivity
contrast(emmeans(mod_old, ~ status_fac * committee),
          list("denied" = c(-1, 1, 0, 0),
               "granted" = c(0, 0, -1, 1)))

data.frame(emmeans(mod_new,  ~ committee_s))[, 2] + fixef(mod_new)$cond[2]

emmeans(mod_old, ~ status_fac * committee)
emmeans(mod_new, ~ status_fac + committee_s + committee)
contrast(emmeans(mod_new, ~ status_fac + committee_s + committee),
         list("denied" = c(1, -1, 0, 0),
              "granted" = c(0, 0, 1, -1)))

options("mesdt.backend" = "glmmTMB")
fit_sdt <- fit_mesdt(formula_mu = ~ emp_gender + (1 | id),
                     formula_lambda = ~ committee + (committee | id),
                     data = dat_exp_2,
                     dv = "assessment",
                     trial_type_var = "status_fac")

unique(dat_exp_1$committee)


#ref_grid_obj <- ref_grid(fit_sdt$fit_obj,
#                         at = list(
#                           committee = unique(dat$committee),
#                           B = c("yes", "no")
#                         ),
#                         cov.reduce = FALSE,
#                         data = mydata,
#                         options = list(contrast = list(A = "contr.treatment", B = "contr.treatment"))
#)



mm_mesdt <- construct_modelmatrices(~ committee, ~ contingencies,
                                    data = dat_exp_2,
                                    trial_type_var = "status_fac",
                                    dv = "response")
mm_mesdt

m <- model.matrix(~ committee + status_fac * contingencies, data = dat_exp_2)


# Marginal Means for Sensitivity
dat_exp_1$status_fac <- dat_exp_1$status_ef
test_m <- glmmTMB(assessment ~ status_ef * emp_gender + (1 | id),
                  data = dat_exp_2,
                  family = binomial("probit"))

emm_sens <- contrast(emmeans(test_m, ~ status_ef * emp_gender),
                     list("male" = c(-1, 1, 0, 0),
                          "female" = c(0, 0, -1, 1)))

emm_sens

summary(test_m)
fixef(test_m)$cond[2] + fixef(test_m)$cond[4]
fixef(test_m)$cond[2] - fixef(test_m)$cond[4]

vcov(test_m)
# Var(X+Y)=Var(X) + Var(Y) + 2Cov(X,Y)
sqrt(vcov(test_m)$cond[2, 2] + vcov(test_m)$cond[4, 4] + 2 * vcov(test_m)$cond[2, 4])
sqrt(vcov(test_m)$cond[2, 2] + vcov(test_m)$cond[4, 4] - 2 * vcov(test_m)$cond[2, 4])

# you don't actually need the other elements of the variance-covariance matrix!
# -> provide to separate emmeans methods for sensitivity and response bias
# using the model matrices etc. given in the model

# example model
options("mlsdt.backend" = "lme4")
fit_sdt <- fit_mesdt(formula_mu = ~ emp_gender + (1 | id),
                     formula_lambda = ~ committee + (committee | id),
                     data = dat_exp_2,
                     dv = "assessment",
                     trial_type_var = "status_fac")

# Try this
recover_data.mesdt <- emmeans:::recover_data.merMod

# rank-deficient models abfangen!
# new argument: sdt_par = c("response_bias", "sensitivity")
# store model matrices in object
# object: whole fit object
# terms: terms of response bias formula
# xlev:
emm_basis.mesdt_fit <- function(object, trms, xlev, grid, ...) {

  if (sdt_par == "sensitivity") {
    m = object$internal$m_frame$mu
    X = object$internal$mm$mu
    # get coefficients for lambda only
    bhat = coef(object)
    # only select rows and columns relevant for sensitivity
    V <- vcov(object$fit_obj)[grep("mu", rownames(vcov(object$fit_obj))),
                                  grep("mu", colnames(vcov(object$fit_obj)))]

    nbasis = matrix(NA)
    # if matrix is rank deficient:
    dfargs <- list(df = summary(object)$d_coef[2])
    dffun = function(k, dfargs) dfargs$df
  } else {
    m = object$internal$m_frame$lambda
    X = object$internal$mm$lambda
    # get coefficients for lambda only
    bhat = coef(object)
    # only select rows and columns relevant for sensitivity
    V <- vcov(object$fit_obj)[grep("lambda", rownames(vcov(object$fit_obj))),
                              grep("lambda", colnames(vcov(object$fit_obj)))]

    nbasis = matrix(NA)
    # if matrix is rank deficient:
    dfargs <- list(df = summary(object)$c_coef[2])
    dffun = function(k, dfargs) dfargs$df
  }

  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs)
}

emmeans(fit_sdt, ~ committee, sdt_par = "response_bias")


#------------------------------------------------------------------------------#
#### Try the same with mesdt ####


options("mesdt.backend" = "glmmTMB")
mod_mesdt <- fit_mesdt(~ committee + (1 | id),
                       ~ committee + (1 | id),
                       dv = "assessment",
                       trial_type_var = "status_fac",
                       data = dat_test)

summary(mod_mesdt$fit_obj)

emmeans(mod_mesdt$fit_obj, ~ 1)



#------------------------------------------------------------------------------#
#### simulate() / predict() ####


test <- simulate(mod_mesdt$fit_obj)
test <- predict(mod_mesdt$fit_obj)

predict(mod_old)


# TODO: test if this works
simulate.mesdt <- function(mesdt_obj, ...) {
  # get method for correct backend
  if (mesdt_obj.backend = "lme4") pred <- lme4::simulate(mesdt_obj$fit_obj, ...)
  else if (mesdt_obj.backend = "glmmTMB") pred <- glmmTMB::simulate(mesdt_obj$fit_obj, ...)

  return(pred)
}


predict.mesdt <- function(mesdt_obj, ...) {
  # get method for correct backend
  if (mesdt_obj.backend = "lme4") pred <- lme4::predict(mesdt_obj$fit_obj, ...)
  else if (mesdt_obj.backend = "glmmTMB") pred <- glmmTMB::predict(mesdt_obj$fit_obj, ...)

  return(pred)
}


