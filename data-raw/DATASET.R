## code to prepare `DATASET` dataset goes here


set.seed(1340569)
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
n_trials <- 200
sigma_mu <- 0.1
sigma_lambda <- 0.1
sigma_mu_diff <- 0.1
sigma_lambda_diff <- 0.1

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
X1 <- c(rep(1, n_subj * n_trials * 2),
        rep(-1, n_subj * n_trials * 2))

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
sim_data$trial_type <- factor(sim_data$trial_type)
sim_data$x1 <- factor(sim_data$x1)
sim_data$ID <- factor(sim_data$ID)

contrasts(sim_data$y) <- contr.sum(2)
contrasts(sim_data$x1) <- contr.sum(2)
contrasts(sim_data$trial_type) <- contr.sum(2)

internal_sdt_data <- sim_data


#------------------------------------------------------------------------------#
#### GLMM for sim data ####

model_test <- glmer(y ~ trial_type * x1 + (trial_type * x1 | ID),
                    family = binomial(link = "probit"), data = sim_data)

model_test_afex <- afex::mixed(y ~ trial_type * x1 + (trial_type * x1 | ID),
                                family = binomial(link = "probit"), data = sim_data,
                                method = "LRT")




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

usethis::use_data(internal_sdt_data, internal_fake_data, model_test, model_test_afex, internal = TRUE, overwrite = T)
#usethis::use_data(DATASET, overwrite = TRUE)
