source("_code_for_all_.R")

set.seed(2024)

# svydesign is available --------------------------------------------------

### check if dr estimators run smoothly --------------------------------------

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
  model_dr_basic_diff <- nonprob(
    selection = ~region  + size,
    outcome = single_shift ~ nace + private,
    svydesign = jvs_svy,
    data = admin,
    pop_size = sum(weights(jvs_svy)),
    method_selection = "logit")
)

expect_silent(
  model_dr_basic_comp <- nonprob(
    selection = ~region  + size,
    outcome = single_shift ~ nace + private,
    svydesign = jvs_svy,
    data = admin,
    pop_size = sum(weights(jvs_svy)),
    control_inference = control_inf(vars_combine = T),
    method_selection = "logit")
)

expect_silent(
  model_dr_basic_boot <- nonprob(
    selection = ~nace,
    outcome = single_shift ~nace,
    svydesign = jvs_svy,
    data = admin,
    control_inference = control_inf(var_method = "bootstrap",
                                    num_boot = 10),
    method_selection = "logit")
)


expect_silent(
  model_dr_basic_boot_comp <- nonprob(
    selection = ~nace + private,
    outcome = single_shift ~nace + size,
    svydesign = jvs_svy,
    data = admin,
    control_inference = control_inf(var_method = "bootstrap",
                                    num_boot = 10,
                                    vars_combine = T),
    method_selection = "logit")
)


expect_silent(
model_dr_varsel <- nonprob(
  selection = ~nace + size,
  outcome = single_shift ~nace + size,
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
  selection = ~nace + size,
  outcome = single_shift ~nace + size,
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
    selection = ~nace,
    outcome = single_shift ~nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T,
                                  vars_combine = T,
                                  bias_correction = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 5)
)
)


# only totals are known ---------------------------------------------------

expect_silent(
  model_dr_basic_totals <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    pop_totals = c(
      "(Intercept)" = sum(weights(jvs_svy)),
      svytotal(~region + private + nace + size, jvs_svy)
    ),
    data = admin,
    method_selection = "logit")
)


expect_silent(
  suppressWarnings(model_dr_basic_totals_boot <- nonprob(
    selection = ~nace + size,
    outcome = single_shift ~nace + size,
    pop_totals = c(
      "(Intercept)" = sum(weights(jvs_svy)),
      svytotal(~nace + size, jvs_svy)
    ),
    data = admin,
    control_inference = control_inf(var_method = "bootstrap",
                                    num_boot = 10),
    method_selection = "logit")
  )
)


## variable combinations
expect_silent(
  suppressWarnings(model_dr_basic_totals_comb <- nonprob(
    selection = ~nace,
    outcome = single_shift ~size,
    pop_totals = c(
      "(Intercept)" = sum(weights(jvs_svy)),
      svytotal(~nace + size, jvs_svy)
    ),
    data = admin,
    control_inference = control_inf(vars_combine = T),
    method_selection = "logit")
  )
)

## variable combinations + boot
expect_silent(
  suppressWarnings(model_dr_basic_totals_comb_boot <- nonprob(
    selection = ~nace,
    outcome = single_shift ~size,
    pop_totals = c(
      "(Intercept)" = sum(weights(jvs_svy)),
      svytotal(~region + private + nace + size, jvs_svy)
    ),
    data = admin,
    control_inference = control_inf(vars_combine = T,
                                    var_method = "bootstrap",
                                    num_boot = 10),
    method_selection = "logit")
  )
)






## TODO FOR LATER: expect different result as IPW weights are not the same as N

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

