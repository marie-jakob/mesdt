test_that("gumbel-min implementation through cloglog gives correct results", {

  library(janitor)
  # TODO: put this into the package
  #dat_exp_2$assessment_rev <- ifelse(dat_exp_2$assessment == "fair", 1, 0)
  #dat_exp_2$status_fac_rev <- ifelse(dat_exp_2$status_fac == 1, -1, 1)

  gumbel_min_mod <- fit_mesdt(~ 1,
                              ~ 1, dv = "assessment", trial_type_var = "status_fac",
                              data = dat_exp_2, distribution = "gumbel-min")
  # Compute this based on aggregated data

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
