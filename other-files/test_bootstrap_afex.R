library(afex)
library(parallel)

# Fit with mesdt
options("mesdt.backend" = "lme4")

fit_mesdt <- fit_mesdt(~ committee + (1 | id),
                       ~ committee + (1 | id),
                 dv = "assessment",
                 trial_type = "status",
                 data = debi3subset)

fit_afex <- mixed(assessment ~ committee * status + (status | id),
                  data = debi3subset,
                  family = binomial("probit"),
                  method = "LRT",
                  test_intercept = FALSE)
logLik(fit_mesdt$fit_obj)
logLik(fit_afex$full_model)

# Same Model

cl <- makeCluster(16, type = "SOCK")
start <- Sys.time()
bootstrap_mesdt <- compute_tests(fit_mesdt, cl = cl,
                                 tests = "bootstrap", nsim = 1000,
                                 test_intercepts = TRUE)
end_mesdt <- Sys.time()
set_backend("glmmTMB")
bootstrap_mesdt_glmmTMB <- compute_tests(fit_mesdt, cl = cl,
                                 tests = "bootstrap", nsim = 1000,
                                 test_intercepts = TRUE)
end_mesdt_tmb <- Sys.time()
start_afex <- Sys.time()
fit_afex <- mixed(assessment ~ committee * status + (status | id),
                  data = debi3subset,
                  family = binomial("probit"),
                  method = "PB",
                  cl = cl)
end_afex <- Sys.time()


