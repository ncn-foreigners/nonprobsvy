# Regression test for the MI bootstrap variance path under method_outcome
# %in% c("nn", "pmm"). Prior to the fix in `boot_mi.R`, the per-replicate
# loop hand-mutated three slots of the `svyrep.design` object and left
# `selfrep` (and `degf`) at their original length. `survey::svytotal()`
# (called from inside method_nn / method_pmm) then failed deterministically
# with the error "(subscript) logical subscript too long". Because the
# loop body sat inside `tryCatch(..., error = function(e) if (verbose)
# message(e))` with no counter increment in the handler, the failure was
# silent and the loop ran forever.
#
# These tests exercise the same code path with a small, fast configuration
# and a global time limit so a regression of the hang would surface as a
# test failure rather than a hung CI run.

suppressPackageStartupMessages(library(survey))

set.seed(2024)
N    <- 1500L
nons <- data.frame(
  x1 = rnorm(300),
  x2 = rnorm(300),
  y  = NA_real_
)
nons$y <- 1 + 0.5 * nons$x1 + 0.3 * nons$x2 + rnorm(300, sd = 0.5)

prob <- data.frame(
  x1 = rnorm(80),
  x2 = rnorm(80),
  w  = rep(N / 80, 80)
)
svy <- survey::svydesign(ids = ~1, weights = ~w, data = prob)

run_boot_quickly <- function(method, num_boot = 5L, time_budget_secs = 60) {
  t0 <- Sys.time()
  res <- nonprob(
    outcome = y ~ x1 + x2,
    data = nons,
    svydesign = svy,
    pop_size = N,
    method_outcome = method,
    family_outcome = "gaussian",
    control_inference = control_inf(
      var_method = "bootstrap",
      num_boot = num_boot,
      cores = 1
    ),
    se = TRUE,
    verbose = FALSE
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  expect_true(elapsed < time_budget_secs,
              info = sprintf("%s MI bootstrap took %.1fs (> %.0fs budget)",
                             method, elapsed, time_budget_secs))
  res
}

# The two formerly broken paths -------------------------------------------
fit_nn  <- run_boot_quickly("nn")
fit_pmm <- run_boot_quickly("pmm")

# Both estimates should be finite scalars and the SE > 0
expect_true(is.finite(fit_nn$output$mean))
expect_true(is.finite(fit_nn$output$SE))
expect_true(fit_nn$output$SE > 0)
expect_true(is.finite(fit_pmm$output$mean))
expect_true(is.finite(fit_pmm$output$SE))
expect_true(fit_pmm$output$SE > 0)

# Confidence interval covers the mean
expect_true(fit_nn$confidence_interval$lower_bound <= fit_nn$output$mean)
expect_true(fit_nn$confidence_interval$upper_bound >= fit_nn$output$mean)
expect_true(fit_pmm$confidence_interval$lower_bound <= fit_pmm$output$mean)
expect_true(fit_pmm$confidence_interval$upper_bound >= fit_pmm$output$mean)

# Multi-outcome MI bootstrap must also run -- exercises the per-outcome
# loop in nonprob_mi.R alongside the bootstrap loop.
nons$y2 <- 2 + 0.3 * nons$x1 - 0.2 * nons$x2 + rnorm(300, sd = 0.5)
t0 <- Sys.time()
fit_multi <- nonprob(
  outcome = y + y2 ~ x1 + x2,
  data = nons,
  svydesign = svy,
  pop_size = N,
  method_outcome = "nn",
  family_outcome = "gaussian",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = 3L, cores = 1),
  se = TRUE,
  verbose = FALSE
)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
expect_true(elapsed < 90,
            info = sprintf("multi-outcome NN bootstrap took %.1fs", elapsed))
expect_true(all(is.finite(fit_multi$output$mean)))
expect_true(all(is.finite(fit_multi$output$SE)))
expect_true(all(fit_multi$output$SE > 0))

# Deterministic per-iteration failure must surface instead of hanging ----
# Directly exercise `boot_mi()` with a stubbed `outcome_method` that always
# errors. With the previous silent tryCatch the loop would never advance;
# the fix surfaces the failure after at most one retry. A 30s budget
# catches a regression of the infinite-loop hang.
ns <- asNamespace("nonprobsvy")
boot_mi_internal     <- get("boot_mi",     envir = ns)
method_nn_original   <- get("method_nn",   envir = ns)

inject_failure <- function(...) stop("forced bootstrap iteration failure")

unlockBinding("method_nn", ns)
assign("method_nn", inject_failure, envir = ns)

mod_obj <- list(
  model          = "nn",
  svydesign      = svy,
  coefficients   = NULL,
  vars_selection = FALSE
)
X_nons_mm <- model.matrix(~ x1 + x2, nons)
X_rand_mm <- model.matrix(~ x1 + x2, prob)

t0 <- Sys.time()
err <- tryCatch(
  boot_mi_internal(
    model_obj        = mod_obj,
    y_nons           = nons$y,
    X_nons           = X_nons_mm,
    X_rand           = X_rand_mm,
    case_weights     = rep(1, nrow(nons)),
    pop_totals       = NULL,
    pop_size         = N,
    family_outcome   = "gaussian",
    control_outcome  = control_out(),
    control_inference= control_inf(var_method = "bootstrap",
                                   num_boot = 3L, cores = 1),
    verbose          = FALSE
  ),
  error = function(e) e
)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

# Restore method_nn no matter what happened
assign("method_nn", method_nn_original, envir = ns)
lockBinding("method_nn", ns)

expect_true(inherits(err, "error"),
            info = "deterministic failure should bubble up as an error")
expect_true(
  grepl("forced bootstrap iteration failure|Bootstrap iteration",
        conditionMessage(err)),
  info = sprintf("unexpected error message: %s", conditionMessage(err))
)
expect_true(elapsed < 30,
            info = sprintf("error surfacing took %.1fs (must not hang)",
                           elapsed))
