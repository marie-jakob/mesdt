library(afex)
library(parallel)

# Fit with mesdt
options("mesdt.backend" = "lme4")
fit_mesdt <- fit_mesdt(~ committee * emp_gender + (1 | id),
                       ~ committee * emp_gender_ef + (1 | id),
                 dv = "assessment",
                 trial_type = "status_fac",
                 data = dat_exp_2)

fit_afex <- mixed(assessment ~ committee * emp_gender * status_fac + (status_fac | id),
                  data = dat_exp_2,
                  family = binomial("probit"),
                  method = "LRT")
logLik(fit_mesdt$fit_obj)
logLik(fit_afex$full_model)

# Same Model

cl <- makeCluster(8, type = "SOCK")
bootstrap_mesdt <- compute_tests(fit_mesdt, cl = cl,
                                 tests = "bootstrap", nsim = 10)

fit_afex <- mixed(assessment ~ committee * emp_gender * status_fac + (status_fac | id),
                  data = dat_exp_2,
                  family = binomial("probit"),
                  method = "PB",
                  cl = cl)



