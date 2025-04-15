source("_code_for_all_.R")

# Create standard IPW estimator
expect_silent(ipw_est1 <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit"
))

# Create calibrated IPW estimator
expect_silent(ipw_est2 <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
))

# population totals

expect_silent(ipw_est1_totals <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit"
))

expect_silent(ipw_est1_means <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_means = pop_means,
  pop_size = N,
  data = admin,
  method_selection = "logit"
))



# Test basic functionality with size variable for standard IPW
expect_silent(result1 <- check_balance(~size, ipw_est1))
expect_silent(result2 <- check_balance(~size, ipw_est2))


# general checks ----------------------------------------------------------

expect_silent(result_dig2 <- check_balance(~size, ipw_est1, dig = 2))
expect_silent(result_dig4 <- check_balance(~size, ipw_est1, dig = 4))

# Should be a list

expect_true(is.list(result1), "Result should be a list")

# Names

expect_true(all(c("nonprob_totals", "prob_totals", "balance") %in% names(result2)),
            "Result should contain all required elements")


# Test error handling - invalid formula
expect_error(check_balance("not_a_formula", ipw_est1),
             "The `x` argument must be a formula")

# Test error handling - missing variable
expect_error(check_balance(~nonexistent_var, ipw_est1),
             "The following variables are not present in the dataset: nonexistent_var.")

# Test error handling - invalid digits
expect_error(check_balance(~size, ipw_est1, dig = -1),
             "The `dig` argument must be a non-negative number")

expect_error(check_balance(~size, ipw_est1, dig = "2"),
             "The `dig` argument must be a non-negative number")

# Test categorical variables (region)
result_cat <- check_balance(~region, ipw_est1)

expect_true(length(result_cat$nonprob_totals) > 1,
            "Categorical variable should have multiple levels")


# Test different digit specifications

expect_true(all(abs(result_dig2$balance - round(result_dig2$balance, 2)) < 1e-10),
            "Rounding to 2 digits should work correctly")
expect_true(all(abs(result_dig4$balance - round(result_dig4$balance, 4)) < 1e-10),
            "Rounding to 4 digits should work correctly")


# check numeric values ----------------------------------------------------

# Test standard IPW values

expect_true(all(c("nonprob_totals", "prob_totals", "balance") %in% names(result1)),
            "Result should contain all required elements")

expect_equal(result1$nonprob_totals,
             c(`(Intercept)` = 52898.1310959942, sizeM = 13529.5497750274, sizeS = 31175.2052120894),
             tolerance = 0.01)

expect_equal(result1$balance,
             c(`(Intercept)` = 1028.13, sizeM = -228.45, sizeS = 1624.21),
             tolerance = 0.01)



# Test calibrated IPW values


expect_equal(result2$nonprob_totals,
             c(`(Intercept)` = 51870.0000000003, sizeM = 13758.0000000002, sizeS = 29550.9999999995),
             tolerance = 0.01)

expect_equal(result2$balance,
             c(`(Intercept)` = 0, sizeM = 0, sizeS = 0),
             tolerance = 0.01)




# test multiple variables ----------------------------------------------------

# Test multiple variables
result_multi <- check_balance(~size + region + nace, ipw_est1)
expect_true(length(result_multi$balance) > 1,
            "Multiple variables should produce multiple balance values")



# Test with multiple categorical variables
result_multi_cat <- check_balance(~region + size + nace, ipw_est1)
expect_true(length(result_multi_cat$balance) > 3,
            "Multiple categorical variables should produce multiple balance values")

# check comparison to provided totals -------------------------------------


expect_silent(
  check_balance(~size, ipw_est1_totals)
)

expect_silent(
  check_balance(~size, ipw_est1_means)
)


expect_equal(
  check_balance(~size, ipw_est1_totals)$balance,
  c(`(Intercept)` = 0, sizeM = 0, sizeS = 0),
  tolerance = 0.01
)
