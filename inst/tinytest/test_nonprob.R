
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
          pop_means = pop_means),
  "The `pop_size` argument must be supplied when `pop_means` is supplied."
)

expect_error(
  nonprob(data = admin,
          outcome = single_shift ~ region,
          target = ~ single_shift,
          pop_size = N)
)

expect_error(
  nonprob(data = admin,
          outcome = single_shift ~ region,
          selection =  ~ region,
          method_outcome = "nn")
)
