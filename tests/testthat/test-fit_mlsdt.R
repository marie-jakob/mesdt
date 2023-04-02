test_that("make_glmer_formula() makes a valid lme4 formula", {
          expect_equal(make_lme_formula(
            form_mu = mu ~ x_test,
            form_lambda = lambda ~ x_test,
            dv = "Y",
            trial_type_var = "trial_type",
            random = "ID",
            within = c("x_test")
            ),
            as.formula("Y ~ trial_type + x_test + x_test:trial_type + (trial_type + x_test + x_test:trial_type | ID)",
                       env = globalenv()))
  }
)

test_that("make_glmer_formula() makes a valid formula without within predictors", {
    expect_equal(make_lme_formula(
      form_mu = mu ~ x_test,
      form_lambda = lambda ~ x_test,
      dv = "Y",
      trial_type_var = "trial_type",
      random = "ID",
      between = c("x_test")
    ),
    as.formula("Y ~ trial_type + x_test + x_test:trial_type + (trial_type | ID)",
              env = globalenv()))
  }
)

