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

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "logit")$output,
  structure(list(mean = 0.72236278190985, SE = 0.042077108910903),
            row.names = "single_shift", class = "data.frame")
)

expect_equal(
  nonprob(
    selection = ~region + private + nace + size,
    target = ~single_shift,
    svydesign = jvs_svy,
    data = admin,
    method_selection = "probit")$output,
  structure(list(mean = 0.723953845426239, SE = 0.0669519861228996),
            row.names = "single_shift", class = "data.frame")
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
  structure(list(mean = 0.722150707253752, SE = 0.042993653604522),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 1e-8
)

expect_equal(
  ipw_cloglog$SE,
  structure(list(prob = 0.0382844151764114, nonprob = 0.0195641970156145),
            row.names = "single_shift", class = "data.frame"),
  tolerance = 1e-8
)
