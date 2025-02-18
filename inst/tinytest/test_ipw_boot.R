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


# standard IPW estimator --------------------------------------------------

# expect_silent(ipw_logit_boot <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "logit",
#   control_inference = control_inf(var_method = "bootstrap", num_boot = 10)
# ))
#
# expect_silent(ipw_logit_boot <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "logit",
#   control_selection = control_sel(est_method = "gee"),
#   control_inference = control_inf(var_method = "bootstrap", num_boot = 10)
# ))

