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

toy_boot_rand <- data.frame(x = c(9, 8), w = c(2, 3))
toy_boot_svy <- svydesign(ids = ~1, weights = ~w, data = toy_boot_rand)

set.seed(11)
fit_nn_exact <- method_nn(
  y_nons = c(0, 50, 100),
  X_nons = matrix(c(0, 5, 10), ncol = 1),
  X_rand = matrix(c(9, 8), ncol = 1),
  svydesign = toy_boot_svy,
  pop_size = 10,
  control_outcome = control_out(k = 1),
  control_inference = control_inf(nn_exact_se = TRUE),
  se = TRUE
)

set.seed(11)
dd_manual <- numeric(50)
for (jj in 1:50) {
  boot_samp <- sample(1:3, size = 3, replace = TRUE)
  y_nons_b <- c(0, 50, 100)[boot_samp]
  X_nons_b <- matrix(c(0, 5, 10)[boot_samp], ncol = 1)
  boot_matches <- RANN::nn2(data = X_nons_b, query = matrix(c(9, 8), ncol = 1), k = 1)
  y_pred_boot <- y_nons_b[boot_matches$nn.idx[, 1]]
  dd_manual[jj] <- stats::weighted.mean(y_pred_boot, c(2, 3))
}

expect_equal(
  fit_nn_exact$var_nonprob,
  var(dd_manual),
  tolerance = 1e-12
)
