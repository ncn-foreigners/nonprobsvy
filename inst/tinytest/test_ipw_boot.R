source("_code_for_all_.R")

# single-core bootstrap ---------------------------------------------------

### standard IPW estimator --------------------------------------------------

## logit
expect_silent(ipw_logit_boot <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
))

expect_equal(dim(ipw_logit_boot$boot_ipw_weights), c(2L, unname(nobs(ipw_logit_boot)["nonprob"])))
expect_true(all(is.finite(ipw_logit_boot$boot_ipw_weights)))
expect_identical(ipw_logit_boot$selection$boot_ipw_weights, ipw_logit_boot$boot_ipw_weights)

expect_silent(suppressWarnings(ipw_logit_totals_boot <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2)
)))

expect_equal(dim(ipw_logit_totals_boot$boot_ipw_weights), c(2L, unname(nobs(ipw_logit_totals_boot)["nonprob"])))
expect_true(all(is.finite(ipw_logit_totals_boot$boot_ipw_weights)))

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
  ipw_logit_boot_multicore <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap", num_boot = 2, cores = 2))
)

expect_equal(dim(ipw_logit_boot_multicore$boot_ipw_weights), c(2L, unname(nobs(ipw_logit_boot_multicore)["nonprob"])))
expect_true(all(is.finite(ipw_logit_boot_multicore$boot_ipw_weights)))

admin_ipw_boot <- admin
admin_ipw_boot$single_shift_inverse <- 1 - admin_ipw_boot$single_shift

expect_silent(
  ipw_mle_totals_boot <- suppressWarnings(nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift + single_shift_inverse,
    pop_totals = pop_totals,
    data = admin_ipw_boot,
    method_selection = "logit",
    control_inference = control_inf(var_method = "bootstrap", num_boot = 2, cores = 2)
  ))
)

expect_equal(dim(ipw_mle_totals_boot$boot_sample), c(2L, 2L))
expect_true(all(is.finite(ipw_mle_totals_boot$output$SE)))
expect_equal(dim(ipw_mle_totals_boot$boot_ipw_weights), c(2L, unname(nobs(ipw_mle_totals_boot)["nonprob"])))
expect_true(all(is.finite(ipw_mle_totals_boot$boot_ipw_weights)))

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

expect_silent(
  ipw_gee_totals_boot <- suppressWarnings(nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee"),
    control_inference = control_inf(var_method = "bootstrap", num_boot = 2, cores = 2)
  ))
)

expect_equal(dim(ipw_gee_totals_boot$boot_sample), c(2L, 1L))
expect_true(all(is.finite(ipw_gee_totals_boot$output$SE)))
