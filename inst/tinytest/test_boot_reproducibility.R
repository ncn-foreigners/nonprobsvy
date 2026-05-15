## Regression tests for issue #115:
## Parallel bootstrap (cores > 1) must be reproducible under set.seed().
## The fix replaced foreach::%dopar% with doRNG::%dorng% in
## R/boot_mi.R, R/boot_ipw.R, and R/boot_dr.R.
##
## CRAN policy caps test cores at 2; running at cores = 2 is sufficient to
## trigger the multicore branch and validate reproducibility.

source("_code_for_all_.R")

n_boot <- 10L
n_cores <- 2L


# 1. IPW bootstrap reproducibility (cores > 1) ----------------------------

set.seed(2024)
ipw_fit1 <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
ipw_fit2 <- nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(ipw_fit1$output$SE, ipw_fit2$output$SE)
expect_identical(ipw_fit1$boot_sample, ipw_fit2$boot_sample)
expect_identical(ipw_fit1$boot_ipw_weights, ipw_fit2$boot_ipw_weights)


# 2. IPW bootstrap reproducibility with pop_totals (cores > 1) ------------

set.seed(2024)
ipw_pt_fit1 <- suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
))

set.seed(2024)
ipw_pt_fit2 <- suppressWarnings(nonprob(
  selection = ~region + private + nace + size,
  target = ~single_shift,
  pop_totals = pop_totals,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
))

expect_identical(ipw_pt_fit1$output$SE, ipw_pt_fit2$output$SE)
expect_identical(ipw_pt_fit1$boot_sample, ipw_pt_fit2$boot_sample)


# 3. MI bootstrap (glm) reproducibility -----------------------------------

set.seed(2024)
mi_glm_fit1 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "glm",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
mi_glm_fit2 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "glm",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(mi_glm_fit1$output$SE, mi_glm_fit2$output$SE)


# 4. MI bootstrap (nn) reproducibility ------------------------------------
# This case exercises NN-specific RNG paths (randomize_returned_ties() etc.)
# inside the workers; before the doRNG fix, this was non-reproducible.

set.seed(2024)
mi_nn_fit1 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "nn",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
mi_nn_fit2 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "nn",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(mi_nn_fit1$output$SE, mi_nn_fit2$output$SE)


# 5. MI bootstrap (pmm) reproducibility -----------------------------------

set.seed(2024)
mi_pmm_fit1 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "pmm",
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
mi_pmm_fit2 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "pmm",
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm"),
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(mi_pmm_fit1$output$SE, mi_pmm_fit2$output$SE)


# 6. Tie-prone NN bootstrap: integer covariates with small k --------------
#   exercises randomize_returned_ties() inside workers.

tie_n <- 60L
tie_x1 <- as.integer(rep(seq(1, 5), length.out = tie_n))
tie_x2 <- as.integer(rep(seq(1, 4), length.out = tie_n))
tie_y  <- tie_x1 + tie_x2 + rnorm(tie_n)
tie_nons <- data.frame(x1 = tie_x1, x2 = tie_x2, y = tie_y)

tie_pop_N <- 5000L
tie_pop_x1 <- as.integer(sample(seq(1, 5), tie_pop_N, replace = TRUE))
tie_pop_x2 <- as.integer(sample(seq(1, 4), tie_pop_N, replace = TRUE))
tie_prob <- data.frame(x1 = tie_pop_x1, x2 = tie_pop_x2, w = rep(50, tie_pop_N))
tie_svy <- svydesign(ids = ~1, weights = ~w, data = tie_prob)

set.seed(2024)
tie_fit1 <- nonprob(
  outcome = y ~ x1 + x2,
  svydesign = tie_svy,
  data = tie_nons,
  method_outcome = "nn",
  control_outcome = control_out(k = 2L),
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
tie_fit2 <- nonprob(
  outcome = y ~ x1 + x2,
  svydesign = tie_svy,
  data = tie_nons,
  method_outcome = "nn",
  control_outcome = control_out(k = 2L),
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(tie_fit1$output$SE, tie_fit2$output$SE)


# 7. DR bootstrap reproducibility -----------------------------------------

set.seed(2024)
dr_fit1 <- nonprob(
  selection = ~region + private + nace + size,
  outcome   = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome   = "glm",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

set.seed(2024)
dr_fit2 <- nonprob(
  selection = ~region + private + nace + size,
  outcome   = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome   = "glm",
  control_inference = control_inf(var_method = "bootstrap",
                                  num_boot = n_boot,
                                  cores = n_cores)
)

expect_identical(dr_fit1$output$SE, dr_fit2$output$SE)
