source("_code_for_all_.R")

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

expect_silent(suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
)))

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

