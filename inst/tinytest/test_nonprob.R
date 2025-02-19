
source("_code_for_all_.R")

# check parameters --------------------------------------------------------

expect_error(
  nonprob(data = admin)
)

expect_error(
  nonprob(data = admin,
          selection = ~ region)
)

expect_error(
  nonprob(data = admin,
          selection = ~ region,
          target = ~ single_shift)
)

expect_error(
  nonprob(data = admin,
          outcome = single_shift ~ region,
          target = ~ single_shift)
)

expect_error(
  nonprob(data = admin,
          outcome = single_shift ~ region,
          target = ~ single_shift,
          pop_means = pop_means)
)

expect_error(
  nonprob(data = admin,
          outcome = single_shift ~ region,
          target = ~ single_shift,
          pop_size = N)
)
