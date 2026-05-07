source("_code_for_all_.R")

toy_rand <- data.frame(x = c(0.5, 4), w = c(2, 3))
toy_svy <- svydesign(ids = ~1, weights = ~w, data = toy_rand)

fit_nn_k1 <- method_nn(
  y_nons = c(1, 4, 9),
  X_nons = matrix(c(0, 2, 5), ncol = 1),
  X_rand = matrix(c(0.5, 4), ncol = 1),
  svydesign = toy_svy,
  pop_size = 10,
  control_outcome = control_out(k = 1),
  se = TRUE
)

expect_equal(
  fit_nn_k1$y_nons_pred,
  c(4, 1, 4)
)

expect_equal(
  fit_nn_k1$var_nonprob,
  301 / 900,
  tolerance = 1e-12
)

expect_equal(
  fit_nn_k1$y_rand_pred,
  c(1, 9)
)

fit_nn_k1_dist <- method_nn(
  y_nons = c(1, 4, 9),
  X_nons = matrix(c(0, 2, 5), ncol = 1),
  X_rand = matrix(c(0.5, 4), ncol = 1),
  svydesign = toy_svy,
  pop_size = 10,
  control_outcome = control_out(k = 1, pmm_weights = "dist"),
  se = TRUE
)

expect_equal(
  fit_nn_k1_dist$y_rand_pred,
  c(1, 9)
)

fit_nonprob_nn_k1 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  method_outcome = "nn",
  data = admin,
  control_outcome = control_out(k = 1)
)

expect_equal(
  fit_nonprob_nn_k1$output,
  data.frame(mean = 0.646597258018859, SE = 0.0279237332302402,
             row.names = "single_shift"),
  tolerance = 1e-6
)
