source("_code_for_all_.R")

set.seed(2024)

expect_stochastic_output <- function(output, expected_mean, mean_tolerance = 1e-12) {
  expect_equal(output$mean, expected_mean, tolerance = mean_tolerance)
  expect_true(is.finite(output$SE) && output$SE > 0)
}

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

set.seed(2024)
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

expect_stochastic_output(model_dr_basic_boot$output, 0.680890499943395)



set.seed(2024)
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

expect_stochastic_output(model_dr_basic_boot_comp$output, 0.704237744834517)

set.seed(2024)
expect_silent(
  model_dr_multi_boot <- nonprob(
    selection = ~nace,
    outcome = dr_y1 + dr_y2 ~ nace,
    svydesign = jvs_svy,
    data = transform(
      admin,
      dr_y1 = as.numeric(single_shift),
      dr_y2 = as.numeric(single_shift) + 0.25 * as.numeric(private)
    ),
    control_inference = control_inf(var_method = "bootstrap", num_boot = 3),
    method_selection = "logit",
    family_outcome = "gaussian"
  )
)

expect_equal(rownames(model_dr_multi_boot$output), c("dr_y1", "dr_y2"))
expect_equal(rownames(model_dr_multi_boot$confidence_interval), c("dr_y1", "dr_y2"))
expect_equal(rownames(model_dr_multi_boot$SE), c("dr_y1", "dr_y2"))
expect_equal(dim(model_dr_multi_boot$boot_sample), c(3L, 2L))
expect_true(all(is.finite(model_dr_multi_boot$output$mean)))
expect_true(all(model_dr_multi_boot$output$SE > 0))
expect_true(all(is.finite(as.matrix(model_dr_multi_boot$confidence_interval))))
expect_true(all(is.na(as.matrix(model_dr_multi_boot$SE))))
expect_equal(
  unname(model_dr_multi_boot$output$SE),
  unname(apply(model_dr_multi_boot$boot_sample, 2, stats::sd)),
  tolerance = 1e-12
)


set.seed(2024)
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
  model_dr_varsel$output$mean,
  0.704444677034764
)

expect_true(is.finite(model_dr_varsel$output$SE) && model_dr_varsel$output$SE > 0)

set.seed(2024)
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

expect_stochastic_output(model_dr_varsel_boot$output, 0.704444677034764)



set.seed(2024)
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


# Kim & Haziza (2014) joint estimator without variable selection ---------------
# Regression test for issue #114: previously crashed with
#   "object 'bias_corr_ys_rand_pred' not found"

set.seed(2024)
expect_silent(
  model_dr_kh <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    pop_size = sum(weights(jvs_svy)),
    method_selection = "logit",
    control_inference = control_inf(bias_correction = TRUE)
  )
)

expect_inherits(model_dr_kh, "nonprob")
expect_true(is.finite(model_dr_kh$output$mean))
expect_true(is.finite(model_dr_kh$output$SE) && model_dr_kh$output$SE > 0)
expect_true(abs(model_dr_kh$output$mean - model_dr_basic$output$mean) > 0)

# bias_correction = TRUE with and without vars_combine should give the same
# joint estimate (the helper is shared between the two code paths).
set.seed(2024)
expect_silent(
  model_dr_kh_comb <- nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin,
    pop_size = sum(weights(jvs_svy)),
    method_selection = "logit",
    control_inference = control_inf(bias_correction = TRUE, vars_combine = TRUE)
  )
)
expect_equal(model_dr_kh$output, model_dr_kh_comb$output)

# Multiple outcomes with bias_correction
set.seed(2024)
expect_silent(
  model_dr_kh_multi <- nonprob(
    selection = ~region + nace + size,
    outcome = dr_y1 + dr_y2 ~ region + nace + size,
    svydesign = jvs_svy,
    data = transform(
      admin,
      dr_y1 = as.numeric(single_shift),
      dr_y2 = as.numeric(single_shift) + 0.25 * as.numeric(private)
    ),
    pop_size = sum(weights(jvs_svy)),
    family_outcome = "gaussian",
    method_selection = "logit",
    control_inference = control_inf(bias_correction = TRUE)
  )
)
expect_equal(rownames(model_dr_kh_multi$output), c("dr_y1", "dr_y2"))
expect_true(all(is.finite(model_dr_kh_multi$output$mean)))
expect_true(all(model_dr_kh_multi$output$SE > 0))

# Bootstrap variance under bias_correction = TRUE
set.seed(2024)
expect_silent(
  model_dr_kh_boot <- nonprob(
    selection = ~region + nace,
    outcome = single_shift ~region + nace,
    svydesign = jvs_svy,
    data = admin,
    pop_size = sum(weights(jvs_svy)),
    method_selection = "logit",
    control_inference = control_inf(bias_correction = TRUE,
                                    var_method = "bootstrap",
                                    num_boot = 10)
  )
)
expect_stochastic_output(model_dr_kh_boot$output, 0.6780598, mean_tolerance = 1e-5)

# bias_correction = TRUE requires svydesign — error when only pop_totals supplied
expect_error(
  nonprob(
    selection = ~region + private + nace + size,
    outcome = single_shift ~region + private + nace + size,
    pop_totals = c(
      "(Intercept)" = sum(weights(jvs_svy)),
      svytotal(~region + private + nace + size, jvs_svy)
    ),
    data = admin,
    method_selection = "logit",
    control_inference = control_inf(bias_correction = TRUE)
  ),
  pattern = "Bias correction.*requires individual-level probability sample data"
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

set.seed(2024)
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

expect_stochastic_output(model_dr_basic_totals_boot$output, 0.704802388170357)


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
set.seed(2024)
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


expect_stochastic_output(model_dr_basic_totals_comb_boot$output, 0.704802388170358)
