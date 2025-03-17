source("_code_for_all_.R")


m1 <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit")


m2 <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  control_inference = control_inf(vars_selection = T),
  control_outcome = control_out(nfolds = 3),
  control_selection = control_sel(nfolds = 3),
  verbose = T
  )

m2a <- nonprob(
  selection = ~region + private + nace + size,
  outcome = single_shift ~region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  family_outcome = "gaussian",
  control_inference = control_inf(vars_selection = T,
                                  bias_correction = T,
                                  vars_combine = T),
  control_outcome = control_out(nfolds = 2),
  control_selection = control_sel(nfolds = 2),
  verbose = T,
  se = FALSE
)

