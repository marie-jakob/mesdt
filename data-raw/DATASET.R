## code to prepare `DATASET` dataset goes here


set.seed(1340569)
library(glmmTMB)

#------------------------------------------------------------------------------#
#### Simulate Test Data ####

# SDT Params:
# Mu_high: 1.2
# Mu_low: 0.8

# lambda_high: -0.5
# lambda_low: -0.1


mu_pop <- 1
mu_diff <- 0.2
lambda_pop <- 0.2
lambda_diff <- 0.3

n_subj <- 20
# number of trials per condition
n_trials <- 50
sigma_mu <- 0.4
sigma_lambda <- 0.4
sigma_mu_diff <- 0.4
sigma_lambda_diff <- 0.4

mu_ind <- rnorm(n_subj, mean = mu_pop, sd = sigma_mu)
mu_diff_ind <- rnorm(n_subj, mean = mu_diff, sd = sigma_mu_diff)

lambda_ind <- rnorm(n_subj, mean = lambda_pop, sd = sigma_lambda)
lambda_diff_ind <- rnorm(n_subj, mean = lambda_diff, sd = sigma_lambda_diff)

mu_ind_high <- mu_ind + mu_diff_ind
mu_ind_low <- mu_ind - mu_diff_ind

lambda_ind_high <- lambda_ind + lambda_diff_ind
lambda_ind_low <- lambda_ind - lambda_diff_ind


signals_high <- numeric(length = n_subj * n_trials)
lures_high <- numeric(length = n_subj * n_trials)

signals_low <- numeric(length = n_subj * n_trials)
lures_low <- numeric(length = n_subj * n_trials)

p_hit_high <- 1 - pnorm(lambda_ind + lambda_diff_ind,
                        mean = mu_ind + mu_diff_ind)
p_hit_low <- 1 - pnorm(lambda_ind - lambda_diff_ind,
                       mean = mu_ind - mu_diff_ind)
p_fa_high <- 1 - pnorm(lambda_ind + lambda_diff_ind, mean = 0)
p_fa_low <-  1 - pnorm(lambda_ind - lambda_diff_ind, mean = 0)


for (i in 1:n_subj) {
  signal_high_tmp <- rbinom(n = n_trials, size = 1, prob = p_hit_high[i])
  signals_high[((i - 1) * n_trials + 1) : (i * n_trials)] <- signal_high_tmp

  signal_low_tmp <- rbinom(n = n_trials, size = 1, prob = p_hit_low[i])
  signals_low[((i - 1) * n_trials + 1) : (i * n_trials)] <- signal_low_tmp


  lure_high_tmp <- rbinom(n = n_trials, size = 1, prob = p_fa_high[i])
  lures_high[((i - 1) * n_trials + 1) : (i * n_trials)] <- lure_high_tmp

  lure_low_tmp <- rbinom(n = n_trials, size = 1, prob = p_fa_low[i])
  lures_low[((i - 1) * n_trials + 1) : (i * n_trials)] <- lure_low_tmp
}

# Data structure:
# 1st half: signals (X1 = 1), 2nd half lures (X1 = -1)

Y <- c(signals_high, signals_low, lures_high, lures_low)
X1 <- c(rep(0.5, n_subj * n_trials * 2),
        rep(-0.5, n_subj * n_trials * 2))

X2 <- c(rep(1, n_subj * n_trials),
        rep(-1, n_subj * n_trials),
        rep(1, n_subj * n_trials),
        rep(-1, n_subj * n_trials))

# IDs for targets
IDs <- c(sort(rep(1:n_subj, n_trials)),
         # IDs for lures
         sort(rep(1:n_subj, n_trials)),
         sort(rep(1:n_subj, n_trials)),
         sort(rep(1:n_subj, n_trials)))

sim3 <- matrix(c(Y, X1, X2, IDs), nrow = n_trials * n_subj * 2 * 2, ncol = 4)

sim_data <- data.frame(sim3)

names(sim_data) <- c("y", "trial_type", "x1", "ID")

sim_data$y <- factor(sim_data$y)
sim_data$trial_type_fac <- sim_data$trial_type * 2
#sim_data$x1 <- factor(sim_data$x1)
sim_data$ID <- factor(sim_data$ID)


internal_sdt_data <- sim_data


#------------------------------------------------------------------------------#
#### GLMM for sim data ####


model_test <- glmer(y ~ x1 * trial_type + (x1 * trial_type | ID),
                    family = binomial(link = "probit"), data = internal_sdt_data, nAGQ = 0)


model_test_afex <- afex::mixed(y ~ x1 * trial_type + (x1 * trial_type | ID),
                                family = binomial(link = "probit"), data = internal_sdt_data,
                                method = "LRT", test_intercept = T)

model_test_uncor <- glmer(y ~ x1 * trial_type + (x1 * trial_type || ID),
                          family = binomial(link = "probit"), data = internal_sdt_data, nAGQ = 0)

model_test_uncor_afex <- afex::mixed(y ~ x1 * trial_type + (x1 * trial_type || ID),
                               family = binomial(link = "probit"), data = internal_sdt_data,
                               method = "LRT", test_intercept = T)

# same with glmmTMB
model_test_tmb <- glmmTMB(y ~ x1 * trial_type + (x1 * trial_type | ID),
                    family = binomial(link = "probit"), data = internal_sdt_data)

model_test_uncor_tmb <- glmmTMB(y ~ x1 * trial_type + (x1 * trial_type | ID),
                    family = binomial(link = "probit"), data = internal_sdt_data)





#------------------------------------------------------------------------------#
#### Data for modeldata stuff ####

d <- data.frame(
  y = factor(c(0, 0, 1, 0, 1, 0, 0, 1, 1, 0)),
  x1 = factor(sample(seq(4, 5), 10, replace = T)),
  x2 = factor(sample(seq(1:3), 10, replace = T)),
  trial_type = factor(sample(0:1, 10, replace = T)),
  ID = factor(seq(1:10))
)


contrasts(d$y) <- contr.sum(2)
contrasts(d$x1) <- contr.sum(2)
contrasts(d$x2) <- contr.sum(3)
contrasts(d$trial_type) <- contr.sum(2)
contrasts(d$ID) <- contr.sum(10)

internal_fake_data <- d


#------------------------------------------------------------------------------#
#### Subset of Data from Exp 2 (Detecting Bias) - LRTs Type II ####


dat_exp_2 <- readRDS("tests/test-dat/dat_exp_2.rds") %>%
  dplyr::filter(id %in% c("regular-4", "regular-5",
                          "reversed-1", "reversed-2",
                          "balanced-1", "balanced-2"))
dat_exp_2 %>%
  dplyr::mutate(committee = factor(ifelse(committee_ef == 1, "granted", "denied")),
                emp_gender = factor(ifelse(emp_gender_ef == 1, "f", "m")),
                # status_fac = factor(status_ef),
                status_fac = status_ef,
                status_ef = status_ef / 2,
                contingencies = factor(contingencies),
                stimulus = sample(c("stim-1", "stim-2", "stim-3",
                                    "stim-4", "stim-5", "stim-6",
                                    "stim-7", "stim-8", "stim-9", "stim-10"),
                                  nrow(dat_exp_2),
                                  replace = T),
                assessment = ifelse(assessment == "fair", 0, 1)) -> dat_exp_2

# contrasts(dat_exp_2$status_fac) <- contr.sum(2)
contrasts(dat_exp_2$committee) <- contr.sum(2)
contrasts(dat_exp_2$emp_gender) <- contr.sum(2)
contrasts(dat_exp_2$contingencies) <- contr.sum(3)

#------------------------------------------------------------------------------#
#### Only Intercepts ####

full_model <- glmer(assessment ~ status_ef + (status_ef | id),
                                 data = dat_exp_2,
                                 family = binomial("probit"),
                                 nAGQ = 0)
fit_red_lambda <- glmer(assessment ~ 0 + status_ef + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)
fit_red_mu <- glmer(assessment ~ 1 + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)

anova(full_model, fit_red_lambda)
anova(full_model, fit_red_mu)

chi_squares_intercepts <- c(
  -2 * (logLik(fit_red_lambda) - logLik(full_model)),
  -2 * (logLik(fit_red_mu) - logLik(full_model))
)


#------------------------------------------------------------------------------#
#### One predictor on mu ####

full_model <- glmer(assessment ~ status_ef + status_ef:committee_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)

# Intercept lambda -> compare with full model
intercept_lambda <- glmer(assessment ~ 0 + status_ef + status_ef:committee_ef + (status_ef | id),
                          data = dat_exp_2,
                          family = binomial("probit"),
                          nAGQ = 0)

# Intercept mu
intercept_mu_full <- glmer(assessment ~ 1 + status_ef + (status_ef | id),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)
intercept_mu_red <- glmer(assessment ~ 1 + (status_ef | id),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)
# committee mu
committee_mu_red <- intercept_mu_full


anova(full_model, intercept_lambda)
anova(intercept_mu_full, intercept_mu_red)
anova(full_model, committee_mu_red)

chi_squares_one_pred_mu_2 <- c(
  2 * (logLik(full_model) - logLik(intercept_lambda)),
  2 * (logLik(intercept_mu_full) - logLik(intercept_mu_red)),
  2 * (logLik(full_model) - logLik(committee_mu_red))
)


#### Type 3
intercept_mu_3 <- glmer(assessment ~ 1 + status_ef:committee_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)

chi_squares_one_pred_mu_3 <- c(
  2 * (logLik(full_model) - logLik(intercept_lambda)),
  2 * (logLik(full_model) - logLik(intercept_mu_3)),
  2 * (logLik(full_model) - logLik(committee_mu_red))
)


#------------------------------------------------------------------------------#
#### One-factorial design ####


full_model <- glmer(assessment ~ status_ef * committee_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)

intercept_lambda_full <- glmer(assessment ~ 1 + status_ef + status_ef:committee_ef + (status_ef | id),
                              data = dat_exp_2,
                              family = binomial("probit"),
                              nAGQ = 0)
intercept_lambda_red <- glmer(assessment ~ 0 + status_ef + status_ef:committee_ef + (status_ef | id),
                               data = dat_exp_2,
                               family = binomial("probit"),
                               nAGQ = 0)
intercept_mu_full <- glmer(assessment ~ committee_ef + status_ef + (status_ef | id),
                               data = dat_exp_2,
                               family = binomial("probit"),
                               nAGQ = 0)
intercept_mu_red <- glmer(assessment ~ committee_ef + (status_ef | id),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)

lambda_committee <- intercept_lambda_full

mu_committee <- intercept_mu_full

anova(intercept_lambda_full, intercept_lambda_red)
anova(intercept_mu_full, intercept_mu_red)
anova(full_model, lambda_committee)
anova(full_model, mu_committee)

chisquares_one_factor_2 <- c(
  -2 * (logLik(intercept_lambda_red) - logLik(intercept_lambda_full)),
  -2 * (logLik(lambda_committee) - logLik(full_model)),
  -2 * (logLik(intercept_mu_red) - logLik(intercept_mu_full)),
  -2 * (logLik(mu_committee) - logLik(full_model))
)


#### Type 3
# -> only different for intercepts

lambda_intercept_3 <-  glmer(assessment ~ 0 + status_ef * committee_ef + (status_ef | id),
                             data = dat_exp_2,
                             family = binomial("probit"),
                             nAGQ = 0)
mu_intercept_3 <-  glmer(assessment ~ committee_ef + committee_ef:status_ef + (status_ef | id),
                             data = dat_exp_2,
                             family = binomial("probit"),
                             nAGQ = 0)

chisquares_one_factor_3 <- c(
  -2 * (logLik(lambda_intercept_3) - logLik(full_model)),
  -2 * (logLik(lambda_committee) - logLik(full_model)),
  -2 * (logLik(mu_intercept_3) - logLik(full_model)),
  -2 * (logLik(mu_committee) - logLik(full_model))
)




#------------------------------------------------------------------------------#
#### Two-factorial design ####


# Compute LRTs manually for testing
# full model - main effects lambda
full_model_main_effects <- glmer(assessment ~ emp_gender_ef + committee_ef + status_ef + status_ef:committee_ef +
                                   status_ef:emp_gender_ef + status_ef:emp_gender_ef:committee_ef + (status_ef | id),
                                 data = dat_exp_2,
                                 family = binomial("probit"),
                                 nAGQ = 0)
# lambda - committee_ef
lambda_committee_red <- glmer(assessment ~ status_ef * emp_gender_ef + status_ef:committee_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
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


# Test intercepts
lambda_intercept_full <- glmer(assessment ~ 1 + status_ef + status_ef:committee_ef * status_ef:emp_gender_ef +
                                 status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                               data = dat_exp_2,
                               family = binomial("probit"),
                               nAGQ = 0)
lambda_intercept_red <- glmer(assessment ~ 0 + status_ef + status_ef:committee_ef * status_ef:emp_gender_ef +
                                status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                              data = dat_exp_2,
                              family = binomial("probit"),
                              nAGQ = 0)
mu_intercept_full <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef + (status_ef | id),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)
mu_intercept_red <- glmer(assessment ~ committee_ef * emp_gender_ef + (status_ef | id),
                          data = dat_exp_2,
                          family = binomial("probit"),
                          nAGQ = 0)


anova(lambda_committee_red, full_model_main_effects)
anova(lambda_emp_gender_red, full_model_main_effects)
anova(model_full, lambda_interaction)
anova(mu_committee_red, full_model_main_effects_mu)
anova(mu_emp_gender_red, full_model_main_effects_mu)
anova(model_full, mu_interaction)

anova(lambda_intercept_red, lambda_intercept_full)
anova(lambda_intercept_red, lambda_intercept_full)

chisquares_two_factors_2 <- c(
  -2 * (logLik(lambda_intercept_red) - logLik(lambda_intercept_full)),
  -2 * (logLik(lambda_committee_red) - logLik(full_model_main_effects)),
  -2 * (logLik(lambda_emp_gender_red) - logLik(full_model_main_effects)),
  -2 * (logLik(lambda_interaction) - logLik(model_full)),
  -2 * (logLik(mu_intercept_red) - logLik(mu_intercept_full)),
  -2 * (logLik(mu_committee_red) - logLik(full_model_main_effects_mu)),
  -2 * (logLik(mu_emp_gender_red) - logLik(full_model_main_effects_mu)),
  -2 * (logLik(mu_interaction) - logLik(model_full))
)

# Compute reduced Type III models

lambda_intercept_3 <- glmer(assessment ~ 0 + committee_ef * emp_gender_ef * status_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)
lambda_committee_3 <- glmer(assessment ~ emp_gender_ef * status_ef + committee_ef:status_ef +
                              committee_ef:emp_gender_ef + committee_ef:emp_gender_ef:status_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)
lambda_emp_gender_3 <- glmer(assessment ~ committee_ef * status_ef + committee_ef:status_ef + emp_gender_ef:status_ef +
                               committee_ef:emp_gender_ef + committee_ef:emp_gender_ef:status_ef + (status_ef | id),
                             data = dat_exp_2,
                             family = binomial("probit"),
                             nAGQ = 0)
lambda_interaction_3 <- glmer(assessment ~ committee_ef * status_ef + emp_gender_ef * status_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)

mu_intercept_3 <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef:committee_ef +
                          status_ef:emp_gender_ef + status_ef:committee_ef:emp_gender_ef + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)
mu_committee_3 <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef +
                          status_ef:emp_gender_ef:committee_ef + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)
mu_emp_gender_3 <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef * committee_ef +
                          status_ef:emp_gender_ef:committee_ef + (status_ef | id),
                        data = dat_exp_2,
                        family = binomial("probit"),
                        nAGQ = 0)
mu_interaction_3 <- glmer(assessment ~ committee_ef * emp_gender_ef + status_ef * committee_ef +
                           status_ef * emp_gender_ef + (status_ef | id),
                         data = dat_exp_2,
                         family = binomial("probit"),
                         nAGQ = 0)




chisquares_two_factors_3 <- c(
  2 * (logLik(model_full) - logLik(lambda_intercept_3)),
  2 * (logLik(model_full) - logLik(lambda_committee_3)),
  2 * (logLik(model_full) - logLik(lambda_emp_gender_3)),
  2 * (logLik(model_full) - logLik(lambda_interaction_3)),
  2 * (logLik(model_full) - logLik(mu_intercept_3)),
  2 * (logLik(model_full) - logLik(mu_committee_3)),
  2 * (logLik(model_full) - logLik(mu_emp_gender_3)),
  2 * (logLik(model_full) - logLik(mu_interaction_3))
)


#------------------------------------------------------------------------------#
#### Factors with > 2 levels ####

model_full <- glmer(assessment ~ contingencies_ef_1 * status_ef + + contingencies_ef_2 * status_ef +
                      (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)
# Type II SS
model_lambda_intercept <- glmer(assessment ~ 1 + status_ef + status_ef:contingencies_ef_1 +
                                  status_ef:contingencies_ef_2 + (status_ef | id),
                                data = dat_exp_2,
                                family = binomial("probit"),
                                nAGQ = 0)

model_lambda_intercept_red <- glmer(assessment ~ 0 + status_ef + status_ef:contingencies_ef_1 +
                                      status_ef:contingencies_ef_2 + (status_ef | id),
                                data = dat_exp_2,
                                family = binomial("probit"),
                                nAGQ = 0)

model_lambda_cont_red <- model_lambda_intercept

model_mu_intercept <- glmer(assessment ~ contingencies_ef_1 + contingencies_ef_2 + status_ef + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)
model_mu_intercept_red <- glmer(assessment ~ contingencies_ef_1 + contingencies_ef_2 + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)
model_mu_cont_red <- model_mu_intercept

chisquares_contingencies_2 <- c(
  -2 * (logLik(model_lambda_intercept_red) - logLik(model_lambda_intercept)),
  -2 * (logLik(model_lambda_cont_red) - logLik(model_full)),
  -2 * (logLik(model_mu_intercept_red) - logLik(model_mu_intercept)),
  -2 * (logLik(model_mu_cont_red) - logLik(model_full))
)


# Type III

lambda_intercept_3 <- glmer(assessment ~ 0 + status_ef * contingencies_ef_1 +
                              status_ef * contingencies_ef_2 + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)
mu_intercept_3 <- glmer(assessment ~ contingencies_ef_1 + contingencies_ef_2 +
                          status_ef:contingencies_ef_1 + status_ef:contingencies_ef_2 + (status_ef | id),
                            data = dat_exp_2,
                            family = binomial("probit"),
                            nAGQ = 0)


chisquares_contingencies_3 <- c(
  -2 * (logLik(lambda_intercept_3) - logLik(model_full)),
  -2 * (logLik(model_lambda_cont_red) - logLik(model_full)),
  -2 * (logLik(mu_intercept_3) - logLik(model_full)),
  -2 * (logLik(model_mu_cont_red) - logLik(model_full))
)

#------------------------------------------------------------------------------#
#### Uncorrelated mu and lambda ####

model_uncor_sdt <- glmer(assessment ~ committee_ef * status_ef + (1 + committee_ef | id) + (0 + status_ef + status_ef:committee_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)


#------------------------------------------------------------------------------#
#### Crossed random effects ####

fit_cross_intercept <- glmer(assessment ~ status_ef + (status_ef | id) + (0 + status_ef | file_name),
                             data = dat_exp_2, family = binomial("probit"), nAGQ = 0)

fit_cross_slopes <- glmer(assessment ~ committee_ef * status_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                             data = dat_exp_2, family = binomial("probit"), nAGQ = 0)

# LRTs for crossed random effects
full_model <- glmer(assessment ~ status_ef * committee_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                    data = dat_exp_2,
                    family = binomial("probit"),
                    nAGQ = 0)

intercept_lambda_full <- glmer(assessment ~ 1 + status_ef + status_ef:committee_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                               data = dat_exp_2,
                               family = binomial("probit"),
                               nAGQ = 0)
intercept_lambda_red <- glmer(assessment ~ 0 + status_ef + status_ef:committee_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                              data = dat_exp_2,
                              family = binomial("probit"),
                              nAGQ = 0)
intercept_mu_full <- glmer(assessment ~ committee_ef + status_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                           data = dat_exp_2,
                           family = binomial("probit"),
                           nAGQ = 0)
intercept_mu_red <- glmer(assessment ~ committee_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                          data = dat_exp_2,
                          family = binomial("probit"),
                          nAGQ = 0)

lambda_committee <- intercept_lambda_full

mu_committee <- intercept_mu_full

anova(intercept_lambda_full, intercept_lambda_red)
anova(intercept_mu_full, intercept_mu_red)
anova(full_model, lambda_committee)
anova(full_model, mu_committee)

chisquares_cross_2 <- c(
  -2 * (logLik(intercept_lambda_red) - logLik(intercept_lambda_full)),
  -2 * (logLik(lambda_committee) - logLik(full_model)),
  -2 * (logLik(intercept_mu_red) - logLik(intercept_mu_full)),
  -2 * (logLik(mu_committee) - logLik(full_model))
)


#### Type 3
# -> only different for intercepts

lambda_intercept_3 <-  glmer(assessment ~ 0 + status_ef * committee_ef +(committee_ef + status_ef | id) + (0 + status_ef | file_name),
                             data = dat_exp_2,
                             family = binomial("probit"),
                             nAGQ = 0)
mu_intercept_3 <-  glmer(assessment ~ committee_ef + committee_ef:status_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                         data = dat_exp_2,
                         family = binomial("probit"),
                         nAGQ = 0)

chisquares_cross_3 <- c(
  -2 * (logLik(lambda_intercept_3) - logLik(full_model)),
  -2 * (logLik(lambda_committee) - logLik(full_model)),
  -2 * (logLik(mu_intercept_3) - logLik(full_model)),
  -2 * (logLik(mu_committee) - logLik(full_model))
)



#------------------------------------------------------------------------------#
#### Random effects tests ####

full_model <- glmmTMB(assessment ~ committee_ef * status_ef + (status_ef | id),
                    data = dat_exp_2,
                    family = binomial("probit"))

mod_lambda_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (0 + status_ef | id),
                              data = dat_exp_2,
                              family = binomial("probit"))

mod_mu_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (1 | id),
                          data = dat_exp_2,
                          family = binomial("probit"))
anova(full_model, mod_lambda_intercept)
anova(full_model, mod_mu_intercept)

chi_squares_rdm_intercepts <- c(
  -2 * (logLik(mod_lambda_intercept) - logLik(full_model)),
  -2 * (logLik(mod_mu_intercept) - logLik(full_model))
)

# Without correlations


full_model <- glmmTMB(assessment ~ committee_ef * status_ef + (status_ef || id),
                    data = dat_exp_2,
                    family = binomial("probit"))

anova(full_model, mod_lambda_intercept)
anova(full_model, mod_mu_intercept)

chi_squares_rdm_intercepts_uncor <- c(
  -2 * (logLik(mod_lambda_intercept) - logLik(full_model)),
  -2 * (logLik(mod_mu_intercept) - logLik(full_model))
)



# Crossed random effects


fit_cross_slopes <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef + status_ef | id) + (0 + status_ef | file_name),
                          data = dat_exp_2, family = binomial("probit"))

cross_lambda_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (0 + committee_ef + status_ef | id) + (0 + status_ef | file_name),
                                data = dat_exp_2, family = binomial("probit"))
cross_mu_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef | id) + (0 + status_ef | file_name),
                            data = dat_exp_2, family = binomial("probit"))
cross_lambda_committee <- glmmTMB(assessment ~ committee_ef * status_ef + (status_ef | id) + (0 + status_ef | file_name),
                                data = dat_exp_2, family = binomial("probit"))
cross_mu_fn <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef + status_ef | id),
                                data = dat_exp_2, family = binomial("probit"))


anova(fit_cross_slopes, cross_lambda_intercept)
anova(fit_cross_slopes, cross_mu_intercept)
anova(fit_cross_slopes, cross_lambda_committee)
anova(fit_cross_slopes, cross_mu_fn)

chi_squares_rdm_cross <- c(
  -2 * (logLik(cross_lambda_intercept) - logLik(fit_cross_slopes)),
  -2 * (logLik(cross_lambda_committee) - logLik(fit_cross_slopes)),
  -2 * (logLik(cross_mu_intercept) - logLik(fit_cross_slopes)),
  -2 * (logLik(cross_mu_fn) - logLik(fit_cross_slopes))
)

# without correlations

fit_cross_slopes_uncor <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef + status_ef || id) + (0 + status_ef || file_name),
                          data = dat_exp_2, family = binomial("probit"))

cross_lambda_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (0 + committee_ef + status_ef || id) + (0 + status_ef || file_name),
                                data = dat_exp_2, family = binomial("probit"))
cross_mu_intercept <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef || id) + (0 + status_ef || file_name),
                            data = dat_exp_2, family = binomial("probit"))
cross_lambda_committee <- glmmTMB(assessment ~ committee_ef * status_ef + (status_ef || id) + (0 + status_ef || file_name),
                                data = dat_exp_2, family = binomial("probit"))
cross_mu_fn <- glmmTMB(assessment ~ committee_ef * status_ef + (committee_ef + status_ef || id),
                     data = dat_exp_2, family = binomial("probit"))


anova(fit_cross_slopes_uncor, cross_lambda_intercept)
anova(fit_cross_slopes_uncor, cross_mu_intercept)
anova(fit_cross_slopes_uncor, cross_lambda_committee)
anova(fit_cross_slopes_uncor, cross_mu_fn)

chi_squares_rdm_cross_uncor <- c(
  -2 * (logLik(cross_lambda_intercept) - logLik(fit_cross_slopes_uncor)),
  -2 * (logLik(cross_lambda_committee) - logLik(fit_cross_slopes_uncor)),
  -2 * (logLik(cross_mu_intercept) - logLik(fit_cross_slopes_uncor)),
  -2 * (logLik(cross_mu_fn) - logLik(fit_cross_slopes_uncor))
)



#------------------------------------------------------------------------------#
#### Store internal data ####


usethis::use_data(internal_sdt_data, internal_fake_data, model_test, model_test_afex,
                  model_test_uncor, model_test_uncor_afex, dat_exp_2, model_uncor_sdt,
                  model_test_tmb, model_test_uncor_tmb,
                  chi_squares_intercepts, chi_squares_one_pred_mu_2, chi_squares_one_pred_mu_3,
                  chisquares_one_factor_2, chisquares_one_factor_3,
                  chisquares_two_factors_2, chisquares_two_factors_3,
                  chisquares_contingencies_2, chisquares_contingencies_3,
                  fit_cross_intercept, fit_cross_slopes,
                  chisquares_cross_2, chisquares_cross_3,
                  chi_squares_rdm_intercepts, chi_squares_rdm_intercepts_uncor,
                  chi_squares_rdm_cross, chi_squares_rdm_cross_uncor,
                  internal = TRUE, overwrite = T)
#usethis::use_data(DATASET, overwrite = TRUE)




