source("_code_for_all_.R")


# test functions ----------------------------------------------------------

expect_equal(
  control_sel(),
  list(est_method = "mle", gee_h_fun = 1, optimizer = "maxLik",
       maxlik_method = "NR", optim_method = "BFGS", epsilon = 1e-04,
       maxit = 500, trace = FALSE, penalty = "SCAD", a_SCAD = 3.7,
       a_MCP = 3, lambda = -1, lambda_min = 0.001, nlambda = 50,
       nfolds = 10, print_level = 0, start_type = "zero", nleqslv_method = "Broyden",
       nleqslv_global = "dbldog", nleqslv_xscalm = "fixed", dependence = FALSE,
       key = NULL)
)

expect_equal(
  control_inf(),
  list(var_method = "analytic", rep_type = "subbootstrap", vars_selection = FALSE,
       vars_combine = FALSE, bias_correction = FALSE, num_boot = 500,
       alpha = 0.05, cores = 1, keep_boot = TRUE, nn_exact_se = FALSE)
)

expect_error(
  control_inf(vars_selection = c(TRUE, FALSE)),
  "'vars_selection' must be a logical scalar"
)

expect_error(
  control_inf(vars_combine = 1),
  "'vars_combine' must be a logical scalar"
)

expect_error(
  control_inf(bias_correction = NA),
  "'bias_correction' must be a logical scalar"
)

expect_equal(
  control_out(),
  list(epsilon = 1e-8, maxit = 100, trace = FALSE, k = 5, penalty = "SCAD",
       a_SCAD = 3.7, a_MCP = 3, lambda = -1, lambda_min = 0.001, nlambda = 100,
       nfolds = 10, treetype = "kd", searchtype = "standard", pmm_match_type = 1,
       pmm_weights = "none", pmm_k_choice = "none", pmm_reg_engine = "glm",
       npar_loess = list(surface = "direct", statistics = "approximate",
                         trace.hat = "approximate", cell = 0.2, iterations = 4L,
                         iterTrace = FALSE))
)

expect_error(
  control_out(lambda = c(0.01, 0.02)),
  "'lambda' must be a finite numeric scalar"
)

expect_error(
  control_out(lambda = -0.01),
  "'lambda' must be non-negative or -1"
)


# test that control methods actually work ------------------------------------------
