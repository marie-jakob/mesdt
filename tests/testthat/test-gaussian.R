test_that("a basic equal-variance gaussian model works", {

  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  gaussian_mod <- fit_mesdt(~ 1,
                            ~ 1,
                            dv = "assessment", trial_type = "status",
                            data = debi3, distribution = "gaussian")
  debi3 %>%
    group_by(assessment, status) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = fair + unfair) %>%
    mutate(fair_rel = fair / n,
           unfair_rel =unfair / n) -> dat_agg

  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]
  # taken from Meyer-Grant et al. (2025)
  d <- qnorm(ht) - qnorm(fa)
  c <- -qnorm(fa) - d / 2
  s <- summary(gaussian_mod)

  expect_equal(s$d_coef[1], d)
  expect_equal(s$c_coef[1], c)

})


test_that("a basic equal-variance gaussian model with a predictor works", {

  skip_if_not_installed("janitor")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  library(janitor)
  library(tidyr)
  library(dplyr)

  gaussian_mod <- fit_mesdt(~ 1 + committee,
                            ~ 1 + committee,
                            dv = "assessment", trial_type = "status",
                            data = debi3, distribution = "gaussian")
  em_c <- data.frame(emmeans(gaussian_mod, ~ committee, dpar = "response bias"))
  em_d <- data.frame(emmeans(gaussian_mod, ~ committee, dpar = "sensitivity"))

  debi3 %>%
    group_by(committee, assessment, status) %>%
    summarize(n = n()) %>%
    pivot_wider(names_from = c(assessment), values_from = n) %>%
    clean_names() %>%
    mutate(n = fair + unfair) %>%
    mutate(fair_rel = fair / n,
           unfair_rel =unfair / n) -> dat_agg

  ht <- dat_agg$unfair_rel[2]
  fa <- dat_agg$unfair_rel[1]
  d <- qnorm(ht) - qnorm(fa)
  c <- -qnorm(fa) - d / 2

  expect_equal(d, em_d[1, 2])
  expect_equal(c, em_c[1, 2])

  ht <- dat_agg$unfair_rel[4]
  fa <- dat_agg$unfair_rel[3]
  d <- qnorm(ht) - qnorm(fa)
  c <- -qnorm(fa) - d / 2


  expect_equal(d, em_d[2, 2])
  expect_equal(c, em_c[2, 2])

})
