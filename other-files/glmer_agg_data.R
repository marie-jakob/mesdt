
# fix GLMM on non-aggregated data (long format)

library(tidyverse)
library(glmmTMB)

dat <- readRDS("other-files/data_prep.rds") %>%
  mutate(assessment = factor(assessment))

dat %>%
  #group_by(assessment, status_ef, committee, vp_gender, id) %>%
  group_by(assessment, status_ef, committee_ef, emp_gender_ef, id) %>%
  summarize(sum = n()) %>%
  pivot_wider(names_from = assessment, values_from = sum) %>%
  mutate(weights = fair + unfair,
         unfair_rel = unfair / weights) -> dat_agg


glmm_form <- formula(assessment ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (status_ef | id))
glmm_form_agg <- formula(unfair_rel ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (status_ef | id))
start_1 <- Sys.time()
glmm_non_agg <- glmer(glmm_form,
                      data = dat,
                      family = binomial("probit"))
end_1 <- Sys.time()


start_2 <- Sys.time()
glmm_agg <- glmer(glmm_form_agg,
                  data = dat_agg,
                  family = binomial("probit"),
                  weights = dat_agg$weights)
end_2 <- Sys.time()

start_3 <- Sys.time()
glmm_non_agg <- glmmTMB(glmm_form,
                        data = dat,
                        family = binomial("probit"))
end_3 <- Sys.time()


start_4 <- Sys.time()
glmm_agg <- glmmTMB(glmm_form_agg,
                    data = dat_agg,
                    family = binomial("probit"),
                    weights = dat_agg$weights)
end_4 <- Sys.time()



simulate(glmm_agg)


logLik(glmm_non_agg)
logLik(glmm_agg)
all(fixef(glmm_non_agg)$cond == fixef(glmm_agg)$cond)
all(unlist(VarCorr(glmm_non_agg)) == unlist(VarCorr(glmm_agg)))

# Huge difference for lme4!
# glmer non-aggregated:
# Time difference of 22.36022 secs
# glmer aggregated:
# Time difference of 0.5181339 secs
# glmmTMB non-aggregated:
# Time difference of 2.147163 secs
# end_4 - start_4
# Time difference of 1.362794 secs
# -> Bei lme4 umtransformieren (insgesamt die schnellste Methode)


#------------------------------------------------------------------------------#
#### Crossed random effects ####

# TODO: test this

dat %>%
  #group_by(assessment, status_ef, committee, vp_gender, id) %>%
  group_by(assessment, status_ef, committee_ef, emp_gender_ef, id, file_name) %>%
  mutate(file_name = factor(file_name)) %>%
  summarize(sum = n()) %>%
  pivot_wider(names_from = assessment, values_from = sum, values_fill = 0) %>%
  mutate(weights = fair + unfair,
         unfair_rel = unfair / weights) -> dat_agg


# -> Weights are not defined
glmm_form <- formula(assessment ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (1 | id) + (1 | file_name))
glmm_form_agg <- formula(unfair_rel ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (1 | file_name))
start_1 <- Sys.time()
glmm_non_agg <- glmer(glmm_form,
                      data = dat,
                      family = binomial("probit"))
end_1 <- Sys.time()


start_2 <- Sys.time()
glmm_agg <- glmer(glmm_form_agg,
                  data = dat_agg,
                  family = binomial("probit"),
                  weights = dat_agg$weights)
end_2 <- Sys.time()

start_3 <- Sys.time()
glmm_non_agg <- glmmTMB(glmm_form,
                        data = dat,
                        family = binomial("probit"))
end_3 <- Sys.time()


start_4 <- Sys.time()
glmm_agg <- glmmTMB(glmm_form_agg,
                    data = dat_agg,
                    family = binomial("probit"),
                    weights = dat_agg$weights)
end_4 <- Sys.time()

# Huge difference for lme4!
# glmer non-aggregated:
# Time difference of 26.25029 secs
# glmer aggregated:
# Time difference of 4.778918 secs
# glmmTMB non-aggregated:
# Time difference of 1.14368 mins
# end_4 - start_4
# Time difference of 1.572966 secs
# -> Selbst bei crossed random-effects schneller für lme4 somehow...

