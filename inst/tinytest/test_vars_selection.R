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

expect_equal(ipw_logit$output$mean, 0.66695569885237, tolerance = 0.001)

expect_equal(
  coef(ipw_logit)$coef_sel[,1],
  c(`(Intercept)` = -1.3115892989764, naceD.E = 0.961915550496649,
    naceF = -0.603377050732427, naceG = -0.498229926003555, naceH = -0.684371136226331,
    naceP = 1.27255233616817),
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

expect_equal(ipw_logit_cal$output$mean, 0.6669557, tolerance = 0.001)

expect_equal(
  coef(ipw_logit_cal)$coef_sel[, 1],
  c(`(Intercept)` = -1.31158929897838, naceD.E = 0.961915550449876,
    naceF = -0.603377050714868, naceG = -0.498229926007702, naceH = -0.684371136220024,
    naceP = 1.27255233608576),
  tolerance = 0.001
)

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
