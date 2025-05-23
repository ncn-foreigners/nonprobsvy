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

expect_equal(
  model_dr_basic$output,
  structure(list(mean = 0.703337799202312, SE = 0.0119326877806591),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_diff$output,
  structure(list(mean = 0.710766429083922, SE = 0.0106414704034084),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_comp$output,
  structure(list(mean = 0.703337799202313, SE = 0.0119326877806591),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_boot$output,
  structure(list(mean = 0.680890499943395, SE = 0.00839488525518477),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_boot_comp$output,
  structure(list(mean = 0.704237744834517, SE = 0.0136338580338005),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_varsel$output,
  structure(list(mean = 0.704444677034764, SE = 0.0107959592634666),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_varsel_boot$output,
  structure(list(mean = 0.704444677034764, SE = 0.0336990276388565),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_varsel_bias$output,
  structure(list(mean = 0.705129098423997, SE = 0.010354361334916),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_totals$output,
  structure(list(mean = 0.704179604031582, SE = 0.00464362593778821),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_totals_boot$output,
  structure(list(mean = 0.704802388170357, SE = 0.00375802225031908),
            row.names = "single_shift", class = "data.frame")
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

expect_equal(
  model_dr_basic_totals_comb$output,
  structure(list(mean = 0.704802388170358, SE = 0.00424047207406606),
            row.names = "single_shift", class = "data.frame")
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


expect_equal(
  model_dr_basic_totals_comb_boot$output,
  structure(list(mean = 0.704802388170358, SE = 0.00660950110534792),
            row.names = "single_shift", class = "data.frame")
)


