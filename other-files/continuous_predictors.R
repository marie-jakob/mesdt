
library(tidyverse)
library(janitor)
library(lme4)

debi3 %>%
  group_by(assessment, status, id) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = c(assessment, status), values_from = n, values_fill = 0) %>%
  clean_names() %>%
  mutate(fair_1 = (fair_1 + 0.5) / 161,
         fair_1_2 = (fair_1_2 + 0.5) / 97,
         unfair_1 = (unfair_1 + 0.5) / 161,
         unfair_1_2 = (unfair_1_2 + 0.5) / 97) -> frequencies


frequencies$d_prime <- qnorm(frequencies$unfair_1_2) - qnorm(frequencies$unfair_1)
frequencies$c <- - qnorm(frequencies$unfair_1) - (frequencies$d_prime / 2)

ind_params <- merge(frequencies, debi3 %>% select(id, age) %>% distinct(), by = "id")

cor.test(ind_params$age, ind_params$d_prime)
cor.test(ind_params$age, ind_params$c)


m <- glmer(assessment ~ status * committee * age + (status * committee | id),
           data = debi3, nAGQ = 0, family = binomial("probit"))

m_mesdt <- fit_mesdt(
  ~ committee * age + (committee | id),
  ~ committee * age + (committee | id),
  data = debi3,
  trial_type = "status",
  dv = "assessment"
)
logLik(m)
logLik(m_mesdt$fit_obj)
df.residual(m)
df.residual(m_mesdt$fit_obj)

fixef(m)
summary(m_mesdt)$c_coef
summary(m_mesdt)$d_coef


#------------------------------------------------------------------------------#
#### Simulate data with a continuous predictor ####

frequencies$pred_d <- frequencies$d_prime * 0.5 + rnorm(nrow(frequencies), 0, 0.5)
cor(frequencies$pred_d, frequencies$d_prime)

frequencies$pred_c <- frequencies$c * -0.5 + rnorm(nrow(frequencies), 0, 0.2)
cor(frequencies$pred_c, frequencies$c)
debi3 <- mesdt::debi3
debi3 <- merge(debi3, frequencies %>% select(id, pred_d, pred_c), by = "id")


m <- glmer(assessment ~ status + pred_c + (status | id),
           data = debi3,
           family = binomial("probit"))


m_mesdt <- fit_mesdt(
  ~ pred_d + (1 | id),
  ~ pred_d +  (1 | id),
  data = debi3,
  trial_type = "status",
  dv = "assessment"
)
m_mesdt

emmeans(m_mesdt, ~ pred_d, dpar = "response bias")

compute_tests(m_mesdt, tests_discriminability = ~ pred_d, tests_response_bias = ~ pred_d)

# TODO level 1 predictors
