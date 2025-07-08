library(lme4)
library(glmmTMB)
library(mgcv)


dat_exp_2 <- readRDS("tests/test-dat/dat_exp_2.rds") %>%
  dplyr::mutate(committee = factor(committee_ef),
                emp_gender = factor(emp_gender_ef),
                status = factor(status_ef),
                contingencies = factor(contingencies))

dat_exp_2_small <- dat_exp_2 %>%
  dplyr::filter(id %in% c("regular-4", "regular-5",
                          "reversed-1", "reversed-2",
                          "balanced-1", "balanced-2"))

start <- Sys.time()
fit_lme <- glmer(assessment ~ status_ef * committee_ef + (committee_ef + status_ef || id),
                 data = dat_exp_2,
                 nAGQ = 1,
                 family = binomial("probit"))
stop <- Sys.time()
time_lme <- stop - start
start <- Sys.time()
fit_tmb <- glmmTMB(assessment ~ status_ef * committee_ef + (1 | id) + (0 + status_ef |id) + (0 + committee_ef | id),
                   data = dat_exp_2,
                   family = binomial("probit"))
stop <- Sys.time()
time_tmb <- stop - start
time_lme
time_tmb
# speed-up of 21.5 !
logLik(fit_lme)
logLik(fit_tmb)
# Slightly higher log likelihood for TMB


