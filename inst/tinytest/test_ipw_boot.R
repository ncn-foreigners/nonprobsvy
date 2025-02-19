## testing bootstrap
set.seed(2024)

# Load required data
data(admin)
data(jvs)

# Create objects ----------------------------------------------------------

# Create survey design object
expect_silent(
  jvs_svy <- svydesign(
    ids = ~1,
    weights = ~weight,
    strata = ~size + nace + region,
    data = jvs
  ))

N <- sum(weights(jvs_svy))
pop_totals <- colSums(model.matrix(~region + private + nace + size, jvs)*jvs$weight)
pop_means <- pop_totals[-1]/N


# single-core bootstrap ---------------------------------------------------

### standard IPW estimator --------------------------------------------------

## logit
expect_silent(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
))

expect_silent(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
))

### calibrated IPW estimator --------------------------------------------------

expect_silent(suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"),
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
)))

expect_silent(suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"),
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
)))


# multi-core bootstrap -----------------------------------------------------

### standard IPW estimator --------------------------------------------------

expect_silent(
  nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2, cores = 2))
)

### calibrated IPW estimator --------------------------------------------------

expect_silent(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee"),
    control_inference = control_inf(var_method = "bootstrap", num_boot = 2, cores = 2))
)

