
# fix GLMM on non-aggregated data (long format)

library(tidyverse)
library(glmmTMB)
library(lme4)

dat <- readRDS("other-files/data_prep.rds") %>%
  mutate(assessment = factor(assessment))

dat %>%
  #group_by(assessment, status_ef, committee, vp_gender, id) %>%
  group_by(assessment, status_ef, committee_ef, emp_gender_ef, id) %>%
  summarize(sum = n()) %>%
  pivot_wider(names_from = assessment, values_from = sum, values_fill = 0) %>%
  mutate(weights = fair + unfair,
         unfair_rel = unfair / weights) -> dat_agg


glmm_form <- formula(assessment ~ status_ef * emp_gender_ef + emp_gender_ef * committee_ef + (status_ef | id))
glmm_form_agg <- formula(unfair_rel ~ status_ef * emp_gender_ef + emp_gender_ef * committee_ef + (status_ef | id))

glmm_form <- formula(assessment ~ status_ef * emp_gender_ef + emp_gender_ef * committee_ef +
                       (status_ef | id) + (1 | file_name))

glmm_cross <- glmer(glmm_form,
                    data = dat,
                    family = binomial("probit"))


#glmm_form <- formula(assessment ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (status_ef | id))
#glmm_form_agg <- formula(unfair_rel ~ committee_ef * emp_gender_ef + status_ef * emp_gender_ef + (status_ef | id))


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
glmm_non_agg_tmb <- glmmTMB(glmm_form,
                        data = dat,
                        family = binomial("probit"))
end_3 <- Sys.time()


start_4 <- Sys.time()
glmm_agg_tmb <- glmmTMB(glmm_form_agg,
                    data = dat_agg,
                    family = binomial("probit"),
                    weights = dat_agg$weights)
end_4 <- Sys.time()




#logLik(glmm_non_agg)
#logLik(glmm_agg)

fixef(glmm_non_agg)
fixef(glmm_agg)

fixef(glmm_non_agg_tmb)
fixef(glmm_agg_tmb)


all(unlist(VarCorr(glmm_non_agg)) == unlist(VarCorr(glmm_agg)))

# Huge difference for lme4!

# glmer non-aggregated:
end_1 - start_1 # Time difference of 48.134 secs
# glmer aggregated:
end_2 - start_2 # Time difference of 1.290914 secs
# glmmTMB non-aggregated:
end_3 - start_3 # Time difference of 4.154411 secs
# glmmTMB aggregated:
end_4 - start_4 # Time difference of 2.847024 secs



#------------------------------------------------------------------------------#
#### Crossed random effects ####

dat %>%
  #group_by(assessment, status_ef, committee, vp_gender, id) %>%
  group_by(assessment, status_ef, committee_ef, emp_gender_ef, id, file_name) %>%
  mutate(file_name = factor(file_name)) %>%
  summarize(sum = n()) %>%
  pivot_wider(names_from = assessment, values_from = sum, values_fill = 0) %>%
  mutate(weights = fair + unfair,
         unfair_rel = unfair / weights) -> dat_agg


glmm_form <- formula(assessment ~ status_ef + (1 | id) + (1 | file_name))
glmm_form_agg <- formula(unfair_rel ~ status_ef + (1 | id) + (1 | file_name))


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
glmm_non_agg_tmb <- glmmTMB(glmm_form,
                        data = dat,
                        family = binomial("probit"))
end_3 <- Sys.time()


start_4 <- Sys.time()
glmm_agg_tmb <- glmmTMB(glmm_form_agg,
                    data = dat_agg,
                    family = binomial("probit"),
                    weights = dat_agg$weights)
end_4 <- Sys.time()

# glmer non-aggregated:
end_1 - start_1 # Time difference of 24.16835 secs
# glmer aggregated:
end_2 - start_2 # Time difference of 23.49468 secs
# glmmTMB non-aggregated:
end_3 - start_3 # Time difference of 42.86207 secs
# end_4 - start_4
end_4 - start_4 # Time difference of 42.30016 secs




fixef(glmm_non_agg)
fixef(glmm_agg)

fixef(glmm_non_agg_tmb)
fixef(glmm_agg_tmb)

df.residual(glmm_non_agg)
df.residual(glmm_agg)

df.residual(glmm_non_agg_tmb)
df.residual(glmm_agg_tmb)


# Compare lme4 and glmmTMB for crossed random effects
glmm_form <- formula(assessment ~ status_ef * emp_gender_ef +
                       (status_ef | id) + (status_ef | file_name))


start_1 <- Sys.time()
glmm_non_agg <- glmer(glmm_form,
                      data = dat,
                      family = binomial("probit"))
end_1 <- Sys.time()

start_3 <- Sys.time()
glmm_non_agg_tmb <- glmmTMB(glmm_form,
                            data = dat,
                            family = binomial("probit"))
end_3 <- Sys.time()

end_1 - start_1
end_3 - start_3


#------------------------------------------------------------------------------#
#### Prepare Implementation ####


debi3$status_num <- ifelse(debi3$status == "signal", 1, -1)
form_disc <- ~ emp_gender + (1 | id)
form_resp_bias <- ~ emp_gender + (1 | id)

trial_type <- "status_num"
dv <- "assessment"

debi3$status_num <- ifelse(debi3$status == "signal", 0.5, -0.5)
debi3$dv_num <- ifelse(debi3$assessment == "fair", 0, 1)

mod_non_agg <- fit_mesdt(form_disc,
                         form_resp_bias,
                         trial_type = trial_type,
                         dv = "assessment",
                         data = debi3)
mod_non_agg$internal$glmer_formula
# extract all variables from formulas

# all(aggregate_data(debi3, form_disc, form_resp_bias) == agg_wide)


grouping_vars <- c(unique(all.vars(form_disc), all.vars(form_resp_bias)), "status_num", "dv_num")
by_all <- debi3[, which(names(debi3) %in% grouping_vars )]
test_agg <- aggregate(debi3$assessment, by = by_all, FUN = length)
#test_agg$trial_type_char <- ifelse(test_agg$status_num == 1, "signal", "noise")
#test_agg[[trial_type]] <- NULL
#test_agg$dv_char <- ifelse(test_agg$status_num == 1, "signal", "noise")
#test_agg[[dv]] <- NULL

id_vars <- c(unique(all.vars(form_disc), all.vars(form_resp_bias)), "status_num")

agg_wide <- reshape(
  test_agg,
  idvar = id_vars,          # column(s) identifying unique rows
  timevar = "dv_num",    # column whose values become new columns
  direction = "wide"
)

new_names <- c("x.0", "x.1")

agg_wide[[new_names[1]]][is.na(agg_wide[[new_names[1]]])] <- 0
agg_wide[[new_names[2]]][is.na(agg_wide[[new_names[2]]])] <- 0


agg_wide$weight <- agg_wide[[new_names[1]]] + agg_wide[[new_names[2]]]
agg_wide$dv_agg <- agg_wide[[new_names[2]]] / agg_wide$weight


mod_agg <- glmer(dv_agg ~ status_num * emp_gender + (status_num | id),
                 data = agg_wide,
                 family = binomial("probit"),
                 weights = agg_wide$weight)

fixef(mod_agg)
mod_non_agg


