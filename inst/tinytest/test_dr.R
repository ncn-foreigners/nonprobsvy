source("_code_for_all_.R")

expect_silent(
  nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit")
)

expect_silent(
  nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T, vars_combine = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 10)
  )
)

expect_silent(
  nonprob(
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
  nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T,
                                  vars_combine = T,
                                  var_method = "bootstrap",
                                  num_boot = 2,
                                  bias_correction = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2, nlambda = 5)
)
)

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



