source("_code_for_all_.R")

# unit-level data ---------------------------------------------------------

# standard IPW estimator --------------------------------------------------
# logit
set.seed(2024)

expect_silent(ipw_logit <- nonprob(
  selection = ~nace,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(vars_selection = TRUE),
  control_selection = control_sel(nfolds = 2, nlambda = 5, penalty = "SCAD")
))

expect_equal(ipw_logit$output$mean, 0.667745127491327, tolerance = 0.001)

expect_equal(
  coef(ipw_logit)$coef_sel[,1],
  c(`(Intercept)` = -1.28526068896553, naceF = -0.629705660477215,
    naceG = -0.524558576077641, naceH = -0.71069978392529,
    naceP = 1.24622365453205),
  tolerance = 0.001
)
#
# # probit
# expect_silent(ipw_probit <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "probit",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(ipw_probit$output$mean, 0.696841025257586, tolerance = 0.001)
#
# expect_equal(
#   ipw_probit$selection$coefficients,
#   c(`(Intercept)` = -0.620622264509895, region04 = 0.571567049756265,
#     region14 = -0.41726958304003, sizeS = -0.465530099350961),
#   tolerance = 0.001
# )
#
#
# # cloglog
# expect_silent(ipw_cloglog <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "cloglog",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(ipw_cloglog$output$mean, 0.6867709, tolerance = 0.001)
#
# expect_equal(
#   ipw_cloglog$selection$coefficients,
#   c(`(Intercept)` = -1.11006067151316, region04 = 0.71062700197917,
#     region14 = -0.708681877918808, region16 = 0.658847353295409,
#     region18 = 1.05074447881999, region30 = -0.437736297317337, private = 0.0565357728536718,
#     naceF = -0.452218227200938, naceG = -0.350100417915597, naceH = -0.654228397426253,
#     naceP = 0.994511676088436, sizeS = -0.653143554084672),
#   tolerance = 0.001
# )
#
#
# calibrated IPW ----------------------------------------------------------
# logit

set.seed(2024)

expect_silent(suppressWarnings(ipw_logit_cal <- nonprob(
  selection = ~nace,
  target = ~single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = control_inf(vars_selection = TRUE),
  control_selection = control_sel(nfolds = 2, nlambda = 5,
                                  est_method = "gee", penalty = "SCAD")
)))

expect_equal(ipw_logit_cal$output$mean, 0.667745127491287, tolerance = 0.001)

expect_equal(
  coef(ipw_logit_cal)$coef_sel[, 1],
  c(`(Intercept)` = -1.28526068896608, naceF = -0.629705660494188,
    naceG = -0.524558576060509, naceH = -0.710699783918584,
    naceP = 1.24622365453422),
  tolerance = 0.001
)

for (method_selection in c("probit", "cloglog")) {
  expect_silent(suppressWarnings(
    ipw_varsel_smoke <- nonprob(
      selection = ~nace,
      target = ~single_shift,
      svydesign = jvs_svy,
      data = admin,
      method_selection = method_selection,
      control_inference = control_inf(vars_selection = TRUE),
      control_selection = control_sel(
        lambda = 0.05,
        maxit = 50,
        penalty = "lasso"
      ),
      se = FALSE
    )
  ))

  expect_true(is.finite(ipw_varsel_smoke$output$mean))
  expect_equal(ipw_varsel_smoke$estimator, "ipw")
  expect_true("(Intercept)" %in% names(coef(ipw_varsel_smoke)$coef_sel[, 1]))
}

expect_silent(suppressWarnings(
  mi_varsel_lambda <- nonprob(
    outcome = single_shift ~ nace,
    svydesign = jvs_svy,
    data = admin,
    method_outcome = "glm",
    family_outcome = "binomial",
    control_inference = control_inf(vars_selection = TRUE),
    control_outcome = control_out(lambda = 0.05, penalty = "lasso"),
    se = FALSE
  )
))

expect_true(is.finite(mi_varsel_lambda$output$mean))
expect_equal(mi_varsel_lambda$estimator, "mi")
expect_true(inherits(mi_varsel_lambda$outcome[[1]], "ncvreg"))
expect_false(inherits(mi_varsel_lambda$outcome[[1]], "cv.ncvreg"))
expect_equal(
  mi_varsel_lambda$outcome[[1]]$lambda[2],
  mi_varsel_lambda$control$control_outcome$lambda
)

# Regression test (L3): variable selection on the totals-only path (no
# `svydesign`) must (a) emit the "experimental" warning and (b) run to a finite
# estimate even when cross-validation drops covariates. Previously, once CV
# reduced `pop_totals` to the selected subset, re-extracting the response in
# nonprob_ipw() routed the still-full outcome formula (`y ~ x1 + x2 + x3 + x4`)
# through model_frame_pop(), whose name check failed with "Selection and
# population totals have different names." The admin/jvs example masked this
# because every `nace` level was always selected; simulated data with explicit
# noise covariates (x3, x4) forces a drop and exercises the crash path.
# This also covers the has_pop_totals branch of loss_theta() in
# src/nonprobCV_cpp.cpp. An exact mean is deliberately NOT asserted: penalty
# selection on this totals-only path is experimental and not trustworthy enough
# to pin a reference value.
set.seed(2024)
sim_N_pop <- 2000
sim_x1 <- rnorm(sim_N_pop)
sim_x2 <- rnorm(sim_N_pop)
sim_x3 <- rnorm(sim_N_pop) # noise: unrelated to selection and outcome
sim_x4 <- rnorm(sim_N_pop) # noise
sim_sel_p <- plogis(-1 + 1.5 * sim_x1 + 1.5 * sim_x2)
sim_R <- rbinom(sim_N_pop, 1, sim_sel_p)
sim_y <- 2 + sim_x1 + sim_x2 + rnorm(sim_N_pop)
sim_pop <- data.frame(x1 = sim_x1, x2 = sim_x2, x3 = sim_x3, x4 = sim_x4, y = sim_y)
sim_sample_A <- sim_pop[sim_R == 1, , drop = FALSE]
sim_pop_totals <- c("(Intercept)" = sim_N_pop,
                    x1 = sum(sim_x1), x2 = sum(sim_x2),
                    x3 = sum(sim_x3), x4 = sum(sim_x4))

sim_warnings <- character(0)
ipw_totals_varsel <- withCallingHandlers(
  nonprob(
    selection = ~ x1 + x2 + x3 + x4,
    target = ~ y,
    pop_totals = sim_pop_totals,
    data = sim_sample_A,
    method_selection = "logit",
    control_inference = control_inf(vars_selection = TRUE),
    control_selection = control_sel(nfolds = 2, nlambda = 5, penalty = "SCAD"),
    se = FALSE
  ),
  warning = function(w) {
    sim_warnings <<- c(sim_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
# (a) the experimental totals-path warning fired
expect_true(any(grepl("experimental", sim_warnings)))
# (b) no crash: a finite point estimate was produced
expect_true(is.finite(ipw_totals_varsel$output$mean))

# # probit
# expect_silent(suppressWarnings(
#   ipw_probit_cal <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "probit",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5, est_method = "gee")
# )))
#
# expect_equal(ipw_probit_cal$output$mean, 0.6871912, tolerance = 0.001)
#
# expect_equal(
#   ipw_probit_cal$selection$coefficients,
#   c(`(Intercept)` = -0.618762005698136, region04 = 0.647020513903461,
#     region14 = -0.412549547150702, sizeS = -0.475326850316793),
#   tolerance = 0.001
# )
#
#
# # cloglog
# expect_silent(suppressWarnings(ipw_cloglog_cal <- nonprob(
#   selection = ~region + private + nace + size,
#   target = ~single_shift,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "cloglog",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5, est_method = "gee")
# )))
#
# expect_equal(ipw_cloglog_cal$output$mean, 0.686266, tolerance = 0.001)
#
# expect_equal(
#   ipw_cloglog_cal$selection$coefficients,
#   c(`(Intercept)` = -1.03654771826791, region04 = 0.950040552282991,
#     region14 = -0.719591806146281, region16 = 0.808903101620628,
#     region18 = 0.966483329826716, region30 = -0.419085449473154,
#     private = -0.0576229233965298, naceF = -0.431463189577142, naceG = -0.332709225679068,
#     naceH = -0.575682674584135, naceP = 0.854664925365864, sizeS = -0.630124936675463
#   ),
#   tolerance = 0.001
# )
#
# # DR estimator (with standard IPW) ------------------------------------------------------------
#
# ### logit
# expect_silent(dr_logit <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "logit",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_logit$output$mean, 0.7032599, tolerance = 0.001)
#
# ### probit
#
# expect_silent(dr_probit <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "probit",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_probit$output$mean, 0.7032552, tolerance = 0.001)
#
# ### cloglog
#
# expect_silent(dr_cloglog <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "cloglog",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_cloglog$output$mean, 0.7033419, tolerance = 0.001)
#
# # DR estimator (with calibrared IPW) ------------------------------------------------------------
#
# ### logit
# expect_silent(dr_logit_gee <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "logit",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5, est_method = "gee"),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_logit_gee$output$mean, 0.7039189, tolerance = 0.001)
#
# ### probit
#
# expect_silent(dr_probit_gee <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "probit",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5, est_method = "gee"),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_probit_gee$output$mean, 0.7041368, tolerance = 0.001)
#
# ### cloglog
#
# expect_silent(dr_cloglog_gee <- nonprob(
#   selection = ~region + private + nace + size,
#   outcome = single_shift ~ region + private + nace + size,
#   svydesign = jvs_svy,
#   data = admin,
#   method_selection = "cloglog",
#   method_outcome = "glm",
#   family_outcome = "binomial",
#   control_inference = control_inf(vars_selection = TRUE),
#   control_selection = control_sel(nfolds = 2, nlambda = 5, est_method = "gee"),
#   control_outcome = control_out(nfolds = 2, nlambda = 5)
# ))
#
# expect_equal(dr_cloglog_gee$output$mean, 0.7038575, tolerance = 0.001)
#
#
# # pop data only -----------------------------------------------------------
#
#
