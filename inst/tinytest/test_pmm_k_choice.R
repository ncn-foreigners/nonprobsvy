source("_code_for_all_.R")

fit_pmm_min_var <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  method_outcome = "pmm",
  data = admin,
  control_outcome = control_out(pmm_k_choice = "min_var", pmm_reg_engine = "glm")
)

fit_pmm_k9 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  method_outcome = "pmm",
  data = admin,
  control_outcome = control_out(k = 9, pmm_reg_engine = "glm")
)

fit_pmm_k10 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  method_outcome = "pmm",
  data = admin,
  control_outcome = control_out(k = 10, pmm_reg_engine = "glm")
)

expect_equal(
  fit_pmm_min_var$control$control_outcome$k,
  9
)

expect_equal(
  fit_pmm_min_var$output,
  fit_pmm_k9$output
)

expect_true(
  fit_pmm_min_var$output$SE < fit_pmm_k10$output$SE
)
