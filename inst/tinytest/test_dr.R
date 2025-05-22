source("_code_for_all_.R")

set.seed(2024)

# check if dr estimators run smoothly --------------------------------------

expect_silent(
model_dr_basic <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit")
)

expect_silent(
  model_dr_basic_boot <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    control_inference = control_inf(var_method = "bootstrap",
                                    num_boot = 10),
    method_selection = "logit")
)

### comment out before submission
expect_silent(
  model_dr_basic_boot_multi <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    control_inference = control_inf(var_method = "bootstrap",
                                    num_boot = 8,
                                    cores = 4),
    method_selection = "logit")
)
### comment out before submission

expect_silent(
model_dr_varsel <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T, vars_combine = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 5)
  )
)

expect_silent(
model_dr_varsel_boot <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T,
                                  vars_combine = T,
                                  var_method = "bootstrap",
                                  num_boot = 2),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 5)
)
)

expect_silent(
  model_dr_varsel_bias <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T,
                                  vars_combine = T,
                                  num_boot = 2,
                                  bias_correction = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 5)
)
)


# check values ------------------------------------------------------------
expect_equal(
  model_dr_basic$output,
  structure(list(mean = 0.703337799202312, SE = 0.0119326877806591),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 0.0001
)

expect_equal(
  model_dr_varsel$output,
  structure(list(mean = 0.703352077362419, SE = 0.0107866849078447),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 0.0001
)

expect_equal(
  model_dr_varsel_boot$output,
  structure(list(mean = 0.703483526995628, SE = 0.000513183400724996),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 0.001
)

expect_equal(
  model_dr_varsel_bias$output,
  structure(list(mean = 0.703920481275204, SE = 0.0101666123541262),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 0.0001
)

## expect different result as IPW weights are not the same as N

# expect_silent(
#   m2 <- nonprob(
#     selection = ~region + private + nace + size,
#     outcome = single_shift ~region + private + nace + size,
#     svydesign = jvs_svy,
#     data = admin,
#     method_selection = "logit")
# )

# m2_bootstrap_bias_corr_multi <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   pop_size = sum(weights(jvs_svy)),
#   method_selection = "logit",
#   control_inference = control_inf(vars_selection = T,
#                                   vars_combine = T,
#                                   var_method = "bootstrap",
#                                   num_boot = 8,
#                                   bias_correction = T,
#                                   cores = 4),
#   control_outcome = control_out(nfolds = 2),
#   control_selection = control_sel(nfolds = 2, nlambda = 5),
#   verbose = T
# )



