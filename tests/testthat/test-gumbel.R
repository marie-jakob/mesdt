test_that("gumbel-min implementation through cloglog gives correct results", {


  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  gumbel_min_mod <- fit_mesdt(~ 1,
                              ~ 1,
                              dv = "assessment", trial_type = "status_fac",
                              data = dat_exp_2, distribution = "gumbel-min")

  dat_exp_2 %>%
    group_by(assessment, status_fac) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = x0 + x1) %>%
    mutate(fair_rel = x0 / n,
           unfair_rel = x1 / n) -> dat_agg
  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]
  # taken from Meyer-Grant et al. (2025)
  g <- -log(-log(ht)) + log(-log(fa))
  kappa <- log(-log(fa)) - g / 2

  expect_equal(g, unname(coef(gumbel_min_mod$fit_obj)[2]))
  expect_equal(kappa, unname(coef(gumbel_min_mod$fit_obj)[1]))
})


test_that("gumbel-min implementation works with response variable with -1 and 1", {
  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  dat_exp_2$dv_fac <- factor(ifelse(dat_exp_2$assessment == 0, -1, 1))

  gumbel_min_mod <- fit_mesdt(~ 1,
                              ~ 1,
                              dv = "dv_fac", trial_type = "status_fac",
                              data = dat_exp_2, distribution = "gumbel-min")
  dat_exp_2 %>%
    group_by(assessment, status_fac) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = x0 + x1) %>%
    mutate(fair_rel = x0 / n,
           unfair_rel = x1 / n) -> dat_agg
  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]

  # taken from Meyer-Grant et al. (2025)
  g <- -log(-log(ht)) + log(-log(fa))
  kappa <- log(-log(fa)) - g / 2

  expect_equal(g, unname(coef(gumbel_min_mod$fit_obj)[2]))
  expect_equal(kappa, unname(coef(gumbel_min_mod$fit_obj)[1]))

})


test_that("gumbel-min implementation works with predictors.", {
  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  gumbel_min_mod <- fit_mesdt(~ 1 + committee,
                              ~ 1 + committee,
                              dv = "dv_fac", trial_type = "status_fac",
                              data = dat_exp_2, distribution = "gumbel-min")

  em_c <- data.frame(emmeans(gumbel_min_mod, ~ committee, dpar = "response bias"))
  em_d <- data.frame(emmeans(gumbel_min_mod, ~ committee, dpar = "sensitivity"))

  dat_exp_2 %>%
    group_by(committee, assessment, status_fac) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = x0 + x1) %>%
    mutate(fair_rel = x0 / n,
           unfair_rel = x1 / n) -> dat_agg

  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]
  # taken from Meyer-Grant et al. (2025)
  g <- -log(-log(ht)) + log(-log(fa))
  kappa <- log(-log(fa)) - g / 2

  expect_equal(g, em_d[1, 2])
  expect_equal(kappa, em_c[1, 2])

  ht <- dat_agg$unfair_rel[4]
  fa <- dat_agg$unfair_rel[3]
  # taken from Meyer-Grant et al. (2025)
  g <- -log(-log(ht)) + log(-log(fa))
  kappa <- log(-log(fa)) - g / 2

  expect_equal(g, em_d[2, 2])
  expect_equal(kappa, em_c[2, 2])
})


#------------------------------------------------------------------------------#
#### Gumbel Max ####

test_that("gumbel-max produces reasonable results.", {
  skip_if_not_installed("ordinal")
  library(ordinal)
  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  gumbel_max_mod <- fit_mesdt(~ 1,
                              ~ 1,
                              dv = "dv_fac", trial_type = "status_fac",
                              data = dat_exp_2, distribution = "gumbel-max")

  dat_exp_2 %>%
    group_by(assessment, status_fac) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = x0 + x1) %>%
    mutate(fair_rel = x0 / n,
           unfair_rel = x1 / n) -> dat_agg

  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]

  cloglog = function(x) log(-log(1-x))

  sens <- cloglog(ht) - cloglog(fa)
  bias <- -cloglog(fa) - (sens / 2)

  expect_equal(unname(coef(gumbel_max_mod$fit_obj)[2]),
               sens)
  expect_equal(unname(coef(gumbel_max_mod$fit_obj)[1]),
               bias * -1)

  kappa_max <- - qgumbel(fa) - g_max / 2

  g <- log(-log(ht)) - log(-log(fa))
  kappa <- - log(-log(fa)) + g / 2

})
