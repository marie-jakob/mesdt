# based on: http://dwoll.de/rexrepos/posts/anovaMixed.html#two-way-repeated-measures-anova-rbf-pq-design

library(afex)
library(lme4)
library(lmerTest)
library(tidyverse)

data("sk2011.1")

#------------------------------------------------------------------------------#
#### Simulation Setup ####


set.seed(123)
P     <- 2               # Xb1
Q     <- 2               # Xb2
R     <- 3               # Xw1
S     <- 3               # Xw2
Njklm <- 20              # obs per cell
Njk   <- Njklm*P*Q       # number of subjects
N     <- Njklm*P*Q*R*S   # number of observations
id    <- gl(Njk,         R*S, N, labels=c(paste("s", 1:Njk, sep="")))
Xb1   <- gl(P,   Njklm*Q*R*S, N, labels=c("CG", "T"))
Xb2   <- gl(Q,   Njklm  *R*S, N, labels=c("f", "m"))
Xw1   <- gl(R,             S, N, labels=c("A", "B", "C"))
Xw2   <- gl(S,   1,           N, labels=c("-", "o", "+"))


mu      <- 100
eB1     <- c(-5, 5)
eB2     <- c(-5, 5)
eW1     <- c(-5, 0, 5)
eW2     <- c(-5, 0, 5)
eB1B2   <- c(-5, 5, 5, -5)
eB1W1   <- c(-5, 5, 2, -2, 3, -3)
eB1W2   <- c(-5, 5, 2, -2, 3, -3)
eB2W1   <- c(-5, 5, 2, -2, 3, -3)
eB2W2   <- c(-5, 5, 2, -2, 3, -3)
eW1W2   <- c(-5, 2, 3, 2, 3, -5, 2, -5, 3)
eB1B2W1 <- c(-5, 5, 5, -5, 2, -2, -2, 2, 3, -3, -3, 3)
eB1B2W2 <- c(-5, 5, 5, -5, 2, -2, -2, 2, 3, -3, -3, 3)
eB1W1W2 <- c(-5, 5, 2, -2, 3, -3, 3, -3, -5, 5, 2, -2, 2, -2, 3, -3, -5, 5)
eB2W1W2 <- c(-5, 5, 2, -2, 3, -3, 3, -3, -5, 5, 2, -2, 2, -2, 3, -3, -5, 5)
# no 3rd-order interaction B1xB2xW1xW2

names(eB1)     <- levels(Xb1)
names(eB2)     <- levels(Xb2)
names(eW1)     <- levels(Xw1)
names(eW2)     <- levels(Xw2)
names(eB1B2)   <- levels(interaction(Xb1, Xb2))
names(eB1W1)   <- levels(interaction(Xb1, Xw1))
names(eB1W2)   <- levels(interaction(Xb1, Xw2))
names(eB2W1)   <- levels(interaction(Xb2, Xw1))
names(eB2W2)   <- levels(interaction(Xb2, Xw2))
names(eW1W2)   <- levels(interaction(Xw1, Xw2))
names(eB1B2W1) <- levels(interaction(Xb1, Xb2, Xw1))
names(eB1B2W2) <- levels(interaction(Xb1, Xb2, Xw2))
names(eB1W1W2) <- levels(interaction(Xb1, Xw1, Xw2))
names(eB2W1W2) <- levels(interaction(Xb2, Xw1, Xw2))


muJKLM <- mu +
  eB1[Xb1] + eB2[Xb2] + eW1[Xw1] + eW2[Xw2] +
  eB1B2[interaction(Xb1, Xb2)] +
  eB1W1[interaction(Xb1, Xw1)] +
  eB1W2[interaction(Xb1, Xw2)] +
  eB2W1[interaction(Xb2, Xw1)] +
  eB2W2[interaction(Xb2, Xw2)] +
  eW1W2[interaction(Xw1, Xw2)] +
  eB1B2W1[interaction(Xb1, Xb2, Xw1)] +
  eB1B2W2[interaction(Xb1, Xb2, Xw2)] +
  eB1W1W2[interaction(Xb1, Xw1, Xw2)] +
  eB2W1W2[interaction(Xb2, Xw1, Xw2)]
muId  <- rep(rnorm(Njk, 0, 3), each=R*S)
mus   <- muJKLM + muId
sigma <- 50

Y  <- round(rnorm(N, mus, sigma), 1)
d2 <- data.frame(id, Xb1, Xb2, Xw1, Xw2, Y)
d1 <- aggregate(Y ~ id + Xw1 + Xb1 + Xb2, data=d2, FUN=mean)

#------------------------------------------------------------------------------#
#### One-way RM ANOVA ####

summary(aov(Y ~ Xw1 + Error(id/Xw1), data=d2))

aov_ez(id = "id", dv = "Y", data = d1, within = "Xw1", anova_table = list(correction = "none"))

fitF_p <- lmer(Y ~ Xw1 + (1|id), data=d1)
summary(fitF_p)
anova(fitF_p)

# SDT
aov_SDT_1 <- glmer(assessment ~ status_ef * committee_ef + (1 | id) + (0 + status_ef | id),
                 data = dat_exp_2,
                 nAGQ = 0,
                 family = binomial("probit"))

aov_SDT_2 <- glmer(assessment ~ status_ef * committee_ef + (status_ef | id),
                   data = dat_exp_2,
                   nAGQ = 0,
                   family = binomial("probit"))

aov_SDT_3 <- glmer(assessment ~ status_ef * committee_ef + (1 | id) + (1 | status_ef:id) + (1 | committee_ef:id),
                   data = dat_exp_2,
                   nAGQ = 0,
                   family = binomial("probit"))


#------------------------------------------------------------------------------#
#### Section ####

aov_SDT_bias <- glmer(assessment ~ committee_ef + (1 | id),
                   data = dat_exp_2,
                   nAGQ = 0,
                   family = binomial("probit"))

aov_SDT_sens <- glmer(assessment ~ 0 + status_ef + status_ef:committee_ef + (0 + status_ef | id),
                      data = dat_exp_2,
                      nAGQ = 0,
                      family = binomial("probit"))





sk2011.1 %>%
  group_by(id, inference) %>%
  summarize(response = mean(response)) -> sk2011_one

summary(aov_ez("id", "response", sk2011_one,
               within = c("inference"),
               anova_table = list(es = "pes", correction = "none")))

lme_one <- lmer(response ~ inference + (1 | id), data = sk2011_one)
anova(lme_one)

# SS equal: 7413.3
# Mean SS equal: 2471.1

summary(lme_one)
# -> in lme4, F values are computed based on vcov() (i.e., the curvature of the
# likelihood function)


#------------------------------------------------------------------------------#
#### Two-way RM ANOVA ####

sk2011.1 %>%
  group_by(id, type, inference) %>%
  summarize(response = mean(response)) %>%
  mutate(inference_ef_1 = ifelse(inference == "MP", 1,
                                 ifelse(inference == "DA", -1, 0)),
         inference_ef_2 = ifelse(inference == "MT", 1,
                                 ifelse(inference == "DA", -1, 0)),
         inference_ef_3 = ifelse(inference == "AC", 1,
                                 ifelse(inference == "DA", -1, 0)),
         type_ef = ifelse(type == "reversed", 1, -1)) -> sk2011_two

contrasts(sk2011_two$inference) <- contr.sum(4)
contrasts(sk2011_two$type) <- contr.sum(2)

two_way_1 <- lmer(response ~ inference_ef_1 + type_ef +
                    inference_ef_1:type_ef +
                    (1 | id) + (1 | inference_ef_1:id) + (1 | type_ef:id), data = sk2011_two %>% filter(inference %in% c("MP", "DA")))



two_way_2 <- lmer(response ~ inference_ef_1 + type_ef +
                    inference_ef_1:type_ef +
                    (1 | id) + (0 + inference_ef_1 | id) +
                    (0 + type_ef | id), data = sk2011_two %>% filter(inference %in% c("MP", "DA")))


anova(two_way)

summary(aov_ez("id", "response", sk2011_two,
               within = c("type", "inference"),
               anova_table = list(es = "pes", correction = "none")))

#aov_ez(id = "id", dv = "Y", data = d2, within = c("Xw1", "Xw2"), anova_table = list(correction = "none"))
summary(two_way)
anova(two_way)


aov_SDT <- glmer(assessment ~ status_ef * committee_ef + (1 | id) + (1 | status_ef:id) +
                   (1 | committee_ef:id) + (1 | committee_ef:status_ef:id),
                 data = dat_exp_2 %>% filter(contingencies == "balanced"),
                 nAGQ = 0,
                 family = binomial("probit"))

aov_SDT_2 <- glmer(assessment ~ status_ef * committee_ef * emp_gender_ef + (1 | id) + (0 + status_ef | id) +
                     (1 | committee_ef:id) + (1 | emp_gender_ef:id) + (0 + status_ef | committee_ef:id) + (0 + status_ef | emp_gender_ef:id),
                 data = dat_exp_2 %>% filter(contingencies == "balanced"),
                 nAGQ = 0,
                 family = binomial("probit"))

