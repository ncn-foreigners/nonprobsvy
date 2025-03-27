source("_code_for_all_.R")

# checking methods --------------------------------------------------------

expect_error(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    pop_totals = pop_totals,
    method_outcome = "nn",
    data = admin)
)

expect_error(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    pop_totals = pop_totals,
    method_outcome = "pmm",
    data = admin)
)

expect_error(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    pop_totals = pop_totals,
    method_outcome = "npar",
    data = admin)
)

