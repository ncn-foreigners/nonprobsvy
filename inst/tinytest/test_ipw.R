source("_code_for_all_.R")

# population data only ----------------------------------------------------------------

#### mle with pop totals -> gee

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit")$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee"))$output
)

#### population totals only

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit",
    start_selection = numeric(NROW(pop_totals)))$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit")$output
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "probit",
    start_selection = numeric(NROW(pop_totals)))$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "probit")$output
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "cloglog",
    start_selection = numeric(NROW(pop_totals)))$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "cloglog")$output,
  tolerance = 0.001
)

## population totals vs population means and pop_size)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit")$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_means = pop_means,
    pop_size = N,
    data = admin,
    method_selection = "logit")$output
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "probit")$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_means = pop_means,
    pop_size = N,
    data = admin,
    method_selection = "probit")$output
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "cloglog")$output,
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_means = pop_means,
    pop_size = N,
    data = admin,
    method_selection = "cloglog")$output,
  tolerance = 0.001
)

# unit-level data ---------------------------------------------------------

### if population totals are presented then they should be equal the gee h=1 method when unit-level data are provided

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_size = pop_totals[1],
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"))$output$mean,

  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_totals = pop_totals,
    data = admin,
    method_selection = "logit")$output$mean
)

expect_warning(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    pop_size = pop_totals[1] + 1,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee"),
    se = FALSE),
  "survey-weight denominator"
)

ipw_gee_warning <- suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_size = pop_totals[1] + 1,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"),
  se = FALSE))

expect_equal(ipw_gee_warning$ipw_estimator, "ht")
expect_equal(ipw_gee_warning$ipw_denominator_source, "survey weights")

ipw_mle_hajek <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  se = FALSE
)

expect_equal(ipw_mle_hajek$ipw_estimator, "hajek")
expect_equal(ipw_mle_hajek$ipw_denominator_source, "estimated IPW weights")
expect_equal(ipw_mle_hajek$ipw_denominator, sum(ipw_mle_hajek$ipw_weights))

ipw_totals_ht <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  se = FALSE
)

expect_equal(ipw_totals_ht$ipw_estimator, "ht")
expect_equal(ipw_totals_ht$ipw_denominator_source, "population totals")
expect_equal(ipw_totals_ht$ipw_denominator, as.numeric(pop_totals[1]))

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit")$output,
  structure(list(mean = 0.708322897640164, SE = 0.00984765550279513),
            row.names = "single_shift", class = "data.frame")
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "probit")$output,
  structure(list(mean = 0.708804223550342, SE = 0.0103134442462989),
            row.names = "single_shift", class = "data.frame")
)

## regression test for the probit + GEE analytic variance bug:
## a scalar `ifelse()` in est_method_ipw.R used to truncate the propensity-score
## derivative vectors (ps_nons_der / est_ps_rand_der) to length 1, corrupting the
## probability-sample variance component (var_cov2) and inflating the SE
## (~0.096 instead of the correct ~0.036). Guard both the derivative lengths and
## the resulting SE.
ipw_probit_gee <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "probit",
  control_selection = control_sel(est_method = "gee"))

# derivatives must stay full-length (n_nonprob / n_prob), not truncated to 1
expect_equal(length(ipw_probit_gee$selection$selection_model$ps_nons_der), nrow(admin))
expect_equal(length(ipw_probit_gee$selection$selection_model$est_ps_rand_der), nrow(jvs))

# probit + GEE SE should match the (never-buggy) logit + GEE SE on these data
expect_equal(
  ipw_probit_gee$output,
  structure(list(mean = 0.704160408297455, SE = 0.036409969848563),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 1e-8
)

## regression test for issue #90: cloglog analytic IPW SE should remain stable
ipw_cloglog <- nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "cloglog")

expect_equal(
  ipw_cloglog$output,
  structure(list(mean = 0.708486797147196, SE = 0.0095719012912863),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 1e-8
)

expect_equal(
  ipw_cloglog$SE,
  structure(list(prob = 0.00825111488581828, nonprob = 0.00485184474928417),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 1e-8
)

expect_equal(ipw_cloglog$ipw_estimator, "hajek")
expect_equal(ipw_cloglog$ipw_denominator_source, "estimated IPW weights")
