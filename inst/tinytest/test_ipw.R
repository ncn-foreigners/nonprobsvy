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
    method_selection = "cloglog")$output
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
