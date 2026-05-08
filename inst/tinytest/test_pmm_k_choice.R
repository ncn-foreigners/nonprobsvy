source("_code_for_all_.R")

min_var_nonprob <- data.frame(
  x = c(-1.8177740, -0.3640923, 0.1124220, 1.5320696, 1.5696762, 1.7618691),
  y = c(0.1292877, 1.7150650, 0.4609162, -1.2650612, -3.6868529, -3.4456620)
)
min_var_prob <- data.frame(
  x = c(-1.90154526, -0.08881612, 0.76282111, 1.18186967),
  w = c(3.275379, 1.649224, 1.954543, 1.694877)
)
min_var_svy <- svydesign(ids = ~1, weights = ~w, data = min_var_prob)

fit_pmm_min_var <- nonprob(
  outcome = y ~ x,
  svydesign = min_var_svy,
  method_outcome = "pmm",
  data = min_var_nonprob,
  control_outcome = control_out(pmm_k_choice = "min_var", pmm_reg_engine = "glm")
)

expect_equal(
  fit_pmm_min_var$control$control_outcome$k,
  5L
)

expect_true(
  is.finite(fit_pmm_min_var$output$SE) && fit_pmm_min_var$output$SE > 0
)

toy_pmm_rand <- data.frame(x = c(1.4, 8.3), w = c(2, 3))
toy_pmm_svy <- svydesign(ids = ~1, weights = ~w, data = toy_pmm_rand)
toy_pmm_y <- c(0, 4, 10, 17, 27, 41)
toy_pmm_X <- cbind(1, c(0, 1.1, 2.7, 4.8, 7.9, 12.5))
toy_pmm_X_rand <- cbind(1, toy_pmm_rand$x)
toy_pmm_weights <- c(1, 2, 3, 5, 8, 13)
colnames(toy_pmm_X) <- colnames(toy_pmm_X_rand) <- c("(Intercept)", "x")

set.seed(27)
fit_pmm_weighted_boot <- method_pmm(
  y_nons = toy_pmm_y,
  X_nons = toy_pmm_X,
  X_rand = toy_pmm_X_rand,
  svydesign = toy_pmm_svy,
  weights = toy_pmm_weights,
  pop_size = 10,
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
  control_inference = control_inf(nn_exact_se = TRUE),
  se = TRUE
)

expect_equal(
  fit_pmm_weighted_boot$y_mi_hat,
  sum(weights(fit_pmm_weighted_boot$svydesign) *
        fit_pmm_weighted_boot$svydesign$variables$y_hat_MI) / 10,
  tolerance = 1e-12
)

set.seed(27)
dd_pmm_manual <- numeric(50)
pmm_control_inference <- control_inf(nn_exact_se = FALSE)
for (jj in 1:50) {
  boot_samp <- sample(seq_along(toy_pmm_y), size = length(toy_pmm_y), replace = TRUE, prob = 1 / toy_pmm_weights)
  y_nons_b <- toy_pmm_y[boot_samp]
  X_nons_b <- toy_pmm_X[boot_samp, , drop = FALSE]
  weights_b <- toy_pmm_weights[boot_samp]

  method_results_boot <- method_glm(
    y_nons = y_nons_b,
    X_nons = X_nons_b,
    X_rand = toy_pmm_X_rand,
    svydesign = toy_pmm_svy,
    weights = weights_b,
    pop_size = 10,
    control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
    control_inference = pmm_control_inference,
    se = FALSE
  )

  pmm_results_boot <- method_nn(
    y_nons = y_nons_b,
    X_nons = method_results_boot$y_nons_pred,
    X_rand = method_results_boot$y_rand_pred,
    svydesign = toy_pmm_svy,
    weights = weights_b,
    pop_size = 10,
    control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
    control_inference = pmm_control_inference,
    se = FALSE
  )

  dd_pmm_manual[jj] <- pmm_results_boot$y_mi_hat
}

expect_equal(
  fit_pmm_weighted_boot$var_nonprob,
  var(dd_pmm_manual),
  tolerance = 1e-12
)
