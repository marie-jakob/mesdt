test_that("make_glmer_formula() makes a valid formula without any predictors", {
    expect_equal(make_glmer_formula(
      form_mu = NULL,
      form_lambda = NULL,
      dv = "Y",
      trial_type_var = "trial_type",
      random = "ID",
    ),
    as.formula("Y ~ trial_type +(trial_type | ID)",
               env = globalenv()))
  }
)

test_that("make_glmer_formula() makes a valid formula without between predictors", {
          expect_equal(make_glmer_formula(
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
    expect_equal(make_glmer_formula(
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


test_that("make_glmer_formula() makes a valid formula with between and within predictors", {
    expect_equal(make_glmer_formula(
      form_mu = mu ~ x_between + x_within,
      form_lambda = lambda ~ x_between + x_within,
      dv = "Y",
      trial_type_var = "trial_type",
      random = "ID",
      within = c("x_within"),
      between = c("x_between")
    ),
    as.formula("Y ~ trial_type + x_between + x_within + x_between:trial_type + x_within:trial_type +
               (trial_type + x_within + x_within:trial_type | ID)",
               env = globalenv()))
  }
)


# TODO:

# empty mu and empty lambda

# different predictors mu and lambda

# differentiate between "+" and "*"


