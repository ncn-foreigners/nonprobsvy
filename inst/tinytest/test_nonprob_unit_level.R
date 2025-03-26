source("_code_for_all_.R")


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



# general tests -----------------------------------------------------------

expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "invalid"
  )
)

expect_error(
  nonprob(
    outcome = single_shift ~ region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    method_outcome = "invalid"
  )
)


expect_error(
  nonprob(
    outcome = single_shift ~ region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    method_outcome = "glm",
    family_outcome = "invalid"
  )
)


## check formulas
expect_error(
  nonprob(
    selection = "not_a_formula",
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin
  ),
  "The `selection` argument must be a formula"
)

expect_error(
  nonprob(
    outcome = "not_a_formula",
    svydesign = jvs_svy,
    data = admin
  ),
  "The `outcome` argument must be a formula"
)


expect_error(
  nonprob(
    data = admin,
    svydesign = jvs_svy
  ),
  "Please provide the `selection` or `outcome` argument"
)

expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin
  ),
  "Please provide the `target` argument as a formula"
)


expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~ region,
    svydesign = jvs_svy,
    data = admin
  ),
  "The following variables in target should not be present in selection"
)


bad_weights <- rep(1, nrow(admin) + 1)
expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    case_weights = bad_weights
  ),
  "Length of the `case_weights` argument must match"
)

expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    case_weights = "not_numeric"
  ),
  "The `case_weights` argument must be a numeric vector"
)


expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = NULL
  ),
  "The `data` argument cannot be an empty data frame."
)




