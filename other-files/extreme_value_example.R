# Use data from Bröder

library(tidyverse)
library(brms)

source("other-files/bin-roc-data.R")

colnames(dbin5point) <- outer(c("miss", "hit", "cr", "fa"), c("cr1", "cr2", "cr3", "cr4", "cr5"),
                              FUN = "paste", sep = "_") %>%
  as.vector()


bdin_all <- dbin5point %>%
  as_tibble()
bdin_all$id <- rownames(dbin5point)
bdin_all <- bdin_all %>%
  separate(id, into = c("exp", "pid"), sep = ", ") %>%
  mutate(pid = str_trim(pid))

bdin_all %>%
  filter(exp == "Broder E3") -> d_broder

bdin_all %>%
  pivot_longer(cols = -c(exp, pid),
               names_to = c("type", "baserate"), names_sep = "_") %>%
  pivot_wider(names_from = "type") %>%
  mutate(Nold = hit + miss,
         Nnew = fa + cr) -> d_broeder_long

d_broeder_long %>%
  rowwise() %>%
  do(tibble(
    exp = .$exp,
    pid = .$pid,
    baserate = .$baserate,
    response_type = rep(c("miss", "hit", "cr", "fa"), times = c(.$miss, .$hit, .$cr, .$fa)),
    old = rep(c(1, 1, 0, 0), times = c(.$miss, .$hit, .$cr, .$fa)),
    response = rep(c(0, 1, 0, 1), times = c(.$miss, .$hit, .$cr, .$fa))
  )) %>%
  select(-response_type) %>%
  mutate(old = ifelse(old == 0, -1, old)) %>%
  ungroup() %>%
  mutate(base_rate_dum = ifelse(baserate == "cr1", 1, 0)) -> d_broeder_bernoulli

options("mesdt.backend" = "glmmTMB")
fit_broeder <- fit_mesdt(discriminability = ~ 1 + (1 | pid),
                 bias = ~ 0 + baserate + (0 + baserate | pid),
                 data = d_broeder_bernoulli,
                 dv = "response",
                 trial_type_var = "old",
                 distribution = "gumbel-min")
summary(fit_broeder)
summary(fit_broeder_brms)
emmeans(fit_broeder, ~baserate, dpar = "response bias")
summary(fit_broeder_brms)



# Fit Henrik
source("../gumbel-reanalysis/gumbelmin_dist-stan.R")
source("../gumbel-reanalysis/gumbelbin-stan.R")

sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "pid")
gumbel_priors <- prior(student_t(3, 1, 2), class = "Intercept") +
  sd_priors

iter <- 3000
warmup <- 1000
gumbel_formula_2 <- brmsformula(
  hit | vint(Nold, fa, Nnew) ~ 1 + (1|p|pid),
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = gumbelbin_family
)



fit_broeder_brms <- brm(
  gumbel_formula_2, data = d_broeder_long,
  stanvars = sv_gumbelbin ,
  prior = gumbel_priors,
  init_r = 0.5,
  iter = iter, warmup = warmup,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(fit_broeder_brms)

emmeans(obj, dpar = "response bias")

summary(fit_broeder)
