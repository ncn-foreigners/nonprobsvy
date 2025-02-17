
# Load test data
data(jvs)
data(admin)

# create objects ----------------------------------------------------------

# Setup survey design
expect_silent(
  jvs_svy <- svydesign(
    ids = ~1,
    weights = ~weight,
    strata = ~size + nace + region,
    data = jvs)
)

# IPW estimator
expect_silent(
  ipw_result <- nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit")
)

# MI Estimator

expect_silent(
  mi_result <- nonprob(
    outcome = single_shift ~ region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

# DR estimator

expect_silent(
  dr_result <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~ region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)



# methods for the IPW -----------------------------------------------------


expect_silent(
  summary(ipw_result)
)

expect_silent(
  nobs(ipw_result)
)

expect_silent(
  pop_size(ipw_result)
)

expect_silent(
  residuals(ipw_result)
)
#
expect_silent(
  residuals(ipw_result, type = "pearson")
)

expect_silent(
  residuals(ipw_result, type = "working")
)

expect_silent(
  residuals(ipw_result, type = "deviance")
)

## we should look into the pearsonSTD residuals
# expect_silent(
#   residuals(ipw_result, "pearsonSTD")
# )

expect_silent(
  cooks.distance(ipw_result)
)

expect_silent(
  hatvalues(ipw_result)
)

expect_silent(
  logLik(ipw_result)
)

expect_silent(
  AIC(ipw_result)
)
#
expect_silent(
  BIC(ipw_result)
)
#
expect_silent(
  confint(ipw_result)
)
#
expect_silent(
  vcov(ipw_result)
)
#
expect_silent(
  deviance(ipw_result)
)


# methods for the MI ------------------------------------------------------


expect_silent(
  summary(mi_result)
)
#
expect_silent(
  nobs(mi_result)
)

expect_silent(
  pop_size(mi_result)
)
#
expect_silent(
  residuals(mi_result)
)

expect_silent(
  residuals(mi_result, type = "pearson")
)
#
expect_silent(
  residuals(mi_result, type = "working")
)

expect_silent(
  residuals(mi_result, type = "deviance")
)

## we should look into the pearsonSTD residuals
# expect_error(
#   residuals(mi_result, "pearsonSTD")
# )

if (isTRUE(tolower(Sys.getenv("TEST_NONPROBSVY_MULTICORE_DEVELOPER")) == "true")) {
  expect_silent(
    cooks.distance(mi_result)
  )
  expect_silent(
    hatvalues(mi_result)
  )
}

expect_silent(
  logLik(mi_result)
)

expect_silent(
  AIC(mi_result)
)
#
expect_silent(
  BIC(mi_result)
)
#
expect_silent(
  confint(mi_result)
)
#
expect_silent(
  vcov(mi_result)
)
#
expect_silent(
  deviance(mi_result)
)


# methods for the MI ------------------------------------------------------

expect_silent(
  summary(dr_result)
)
#
expect_silent(
  nobs(dr_result)
)

expect_silent(
  pop_size(dr_result)
)
#
expect_silent(
  residuals(dr_result)
)

expect_silent(
  residuals(dr_result, type = "pearson")
)
#
expect_silent(
  residuals(dr_result, type = "working")
)

expect_silent(
  residuals(dr_result, type = "deviance")
)

## we should look into the pearsonSTD residuals
# expect_error(
#   residuals(dr_result, "pearsonSTD")
# )

if (isTRUE(tolower(Sys.getenv("TEST_NONPROBSVY_MULTICORE_DEVELOPER")) == "true")) {
  expect_silent(
    cooks.distance(dr_result)
  )

  expect_silent(
    hatvalues(dr_result)
  )
}

expect_silent(
  logLik(dr_result)
)

expect_silent(
  AIC(dr_result)
)

expect_silent(
  BIC(dr_result)
)

expect_silent(
  confint(dr_result)
)

expect_silent(
  vcov(dr_result)
)

expect_silent(
  deviance(dr_result)
)


# methods not implemented -------------------------------------------------

expect_error(
  anova(ipw_result),
  "The `anova` method is not implemented for the `nonprob` class"
)

expect_error(
  plot(ipw_result),
  "We do not provide tools for visual assessment of the results"
)

