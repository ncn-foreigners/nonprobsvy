source("_code_for_all_.R")

known_n_svy <- svydesign(
  ids = ~1,
  weights = ~w,
  data = data.frame(x = c(0.5, 4), w = c(2, 3))
)

fit_nn_known_n <- method_nn(
  y_nons = c(1, 4, 9),
  X_nons = matrix(c(0, 2, 5), ncol = 1),
  X_rand = matrix(c(0.5, 4), ncol = 1),
  svydesign = known_n_svy,
  pop_size = 10,
  control_outcome = control_out(k = 1),
  se = FALSE
)

fit_nn_estimated_n <- method_nn(
  y_nons = c(1, 4, 9),
  X_nons = matrix(c(0, 2, 5), ncol = 1),
  X_rand = matrix(c(0.5, 4), ncol = 1),
  svydesign = known_n_svy,
  control_outcome = control_out(k = 1),
  se = FALSE
)

nn_total <- sum(weights(known_n_svy) * fit_nn_known_n$y_rand_pred)

expect_equal(fit_nn_known_n$y_mi_hat, nn_total / 10)
expect_equal(fit_nn_estimated_n$y_mi_hat, nn_total / sum(weights(known_n_svy)))
expect_true(fit_nn_known_n$y_mi_hat != fit_nn_estimated_n$y_mi_hat)

fit_pmm_known_n <- method_pmm(
  y_nons = c(1, 4, 9),
  X_nons = cbind(1, c(0, 2, 5)),
  X_rand = cbind(1, c(0.5, 4)),
  svydesign = known_n_svy,
  pop_size = 10,
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
  se = FALSE
)

fit_pmm_estimated_n <- method_pmm(
  y_nons = c(1, 4, 9),
  X_nons = cbind(1, c(0, 2, 5)),
  X_rand = cbind(1, c(0.5, 4)),
  svydesign = known_n_svy,
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
  se = FALSE
)

pmm_total <- sum(weights(known_n_svy) * fit_pmm_known_n$svydesign$variables$y_hat_MI)

expect_equal(fit_pmm_known_n$y_mi_hat, pmm_total / 10)
expect_equal(fit_pmm_estimated_n$y_mi_hat, pmm_total / sum(weights(known_n_svy)))
expect_true(fit_pmm_known_n$y_mi_hat != fit_pmm_estimated_n$y_mi_hat)
