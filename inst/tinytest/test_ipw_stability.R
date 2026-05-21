extreme_eta <- c(-1000, -100, -40, 0, 40, 100, 1000)
extreme_x <- matrix(extreme_eta, ncol = 1)
unit_weights <- rep(1, length(extreme_eta))

for (link in c("logit", "probit", "cloglog")) {
  method <- method_ps(link)

  ps <- method$make_link_inv(extreme_eta)
  expect_true(all(is.finite(ps)))
  expect_true(all(ps > 0 & ps < 1))
  expect_true(all(is.finite(1 / ps)))

  expect_true(all(is.finite(method$make_link_inv_der(extreme_eta))))
  expect_true(all(is.finite(method$make_link_inv_rev(extreme_eta))))
  expect_true(all(is.finite(method$make_link_inv_rev_der(extreme_eta))))
  expect_true(all(is.finite(method$make_link_inv_rev_der2(extreme_eta))))

  log_like <- method$make_log_like(
    X_nons = extreme_x,
    X_rand = extreme_x,
    weights = unit_weights,
    weights_rand = unit_weights
  )(1)
  expect_true(is.finite(log_like))
}

logit_method <- method_ps("logit")
logit_log_like <- logit_method$make_log_like(
  X_nons = matrix(0, nrow = 1, ncol = 1),
  X_rand = matrix(1000, nrow = 1, ncol = 1),
  weights = 1,
  weights_rand = 1
)(1)

expect_equal(logit_log_like, -1000, tolerance = 1e-12)
