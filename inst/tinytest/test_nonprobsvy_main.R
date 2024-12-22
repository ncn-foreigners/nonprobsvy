# Definition of data
#library(survey)
# library(dplyr)
# data <- read.csv("test_data.csv")
# data |>
#   mutate(zawod_kod2 = factor(zawod_kod2),
#          jedna_zmiana=as.numeric(jedna_zmiana))
# data$zawod_kod2 <- as.factor(data$zawod_kod2)
# data$jedna_zmiana <- as.numeric(data$jedna_zmiana)

# ## probability sample
# data |>
#   filter(is.na(id_popyt)) |>
#   select(id_jednostki, klasa_pr, sek, woj, zawod_kod2, wolne_miejsca_cbop, jedna_zmiana) -> cbop_df
#
# ## nonprobability sample
# data |>
#   filter(!is.na(id_popyt)) |>
#   select(id_jednostki, klasa_pr, sek, woj, zawod_kod2, waga, wolne_miejsca) -> popyt_df

# write.csv(cbop_df, "cbop_df.csv")
# write.csv(popyt_df, "popyt_df.csv")

# popyt_df <- read.csv("popyt_df.csv")
# cbop_df <- read.csv("cbop_df.csv")
#
#
# popyt_svy <- svydesign(ids=~1, strata = ~klasa_pr+sek+woj, weights = ~waga*wolne_miejsca,
#                        data = popyt_df)
# cbop_df_long <- cbop_df[rep(1:nrow(cbop_df),cbop_df$wolne_miejsca_cbop), ]
#
# # IPW
# expect_silent(
#   test_ipw_1 <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                       target = ~ jedna_zmiana ,
#                       data = cbop_df_long,
#                       svydesign = popyt_svy,
#                       method_selection = "logit")
# )
#
# expect_equivalent(test_ipw_1$output$mean,
#                   0.621989, # TODO,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_ipw_1$confidence_interval[1] < 0.621989) & (0.621989 < test_ipw_1$confidence_interval[2])
# )
#
# expect_silent(
#   test_ipw_2 <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                       target = ~ jedna_zmiana,
#                       data = cbop_df,
#                       svydesign = popyt_svy,
#                       weights = cbop_df$wolne_miejsca_cbop)
# )
#
# expect_equivalent(test_ipw_2$output$mean,
#                   0.621989, # TODO,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_ipw_2$confidence_interval[1] < 0.621989) & (0.621989 < test_ipw_2$confidence_interval[2])
# )

# # SCAD
# expect_silent(
#   test_ipw_1_scad <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                         target = ~ jedna_zmiana ,
#                         data = cbop_df_long,
#                         svydesign = popyt_svy,
#                         method_selection = "logit",
#                         control_inference = controlInf(vars_selection = TRUE),
#                         control_selection = controlSel(penalty = "SCAD"))
# )
#
# expect_equivalent(test_ipw_1_scad$output$mean,
#                   0.6219899, # TODO,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_ipw_1_scad$confidence_interval[1] < 0.6219899) & (0.6219899 < test_ipw_1_scad$confidence_interval[2])
# )
#
# # LASSO
# expect_silent(
#   test_ipw_1_lasso <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                         target = ~ jedna_zmiana ,
#                         data = cbop_df_long,
#                         svydesign = popyt_svy,
#                         method_selection = "logit",
#                         control_inference = controlInf(vars_selection = TRUE),
#                         control_selection = controlSel(penalty = "lasso"))
#   )
#
#   expect_equivalent(test_ipw_1_lasso$output$mean,
#                     0.6191018, # TODO,
#                     tolerance = 0.1)
#
#   expect_true(
#     (test_ipw_1_lasso$confidence_interval[1] < 0.6191018) & (0.6191018 < test_ipw_1$confidence_interval[2])
#   )
#
# # MCP
# expect_silent(
#   test_ipw_1_mcp <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                         target = ~ jedna_zmiana ,
#                         data = cbop_df_long,
#                         svydesign = popyt_svy,
#                         method_selection = "logit",
#                         control_inference = controlInf(vars_selection = TRUE),
#                         control_selection = controlSel(penalty = "MCP"))
# )
#
# expect_equivalent(test_ipw_1_mcp$output$mean,
#                   0.621989, # TODO,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_ipw_1_mcp$confidence_interval[1] < 0.6191018) & (0.6191018 < test_ipw_1$confidence_interval[2])
# )

# DR
# expect_silent(
#   test_dr_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy)
# )
#
# expect_equivalent(test_dr_1$output$mean,
#                   0.6048031,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1$confidence_interval[1] < 0.6048031) &
#     (0.6048031< test_dr_1$confidence_interval[2])
# )
#
# pop_totals <- c(`(Intercept)` = 294294, klasa_prM = 128940, klasa_prS = 71230)
# test_dr_2 <- nonprob(selection = ~ klasa_pr,
#                     outcome = jedna_zmiana ~ klasa_pr,
#                     data = cbop_df,
#                     pop_totals = pop_totals,
#                     control_inference = controlInf(var_method = "analytic")) # TODO warning to connected to algorithm convergence
#
# expect_equivalent(test_dr_2$output$mean,
#                   0.6747754,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_2$confidence_interval[1] < 0.6747754) &
#     (0.6747754 < test_dr_2$confidence_interval[2])
# )
#
# # PMM
# test_dr_1_pmm <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "pmm")
#
# expect_equivalent(test_dr_1_pmm$output$mean,
#                   0.6490007,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_pmm$confidence_interval[1] < 0.6490007) &
#     (0.6490007< test_dr_1_pmm$confidence_interval[2])
# )
# # NN
# test_dr_1_nn <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "nn")
#
# expect_equivalent(test_dr_1_nn$output$mean,
#                   0.6003411,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_nn$confidence_interval[1] < 0.6003411) &
#     (0.6003411< test_dr_1_nn$confidence_interval[2])
# )

# bias minimization
# expect_silent(
#   test_dr_1_bm <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        control_inference = controlInf(bias_correction = TRUE))
# )
#
# expect_equivalent(test_dr_1_bm$output$mean,
#                   0.6064716,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_bm$confidence_interval[1] < 0.6064716) &
#     (0.6064716 < test_dr_1_bm$confidence_interval[2])
# )

# # SCAD
# expect_silent(
#   test_dr_1_scad <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_selection = controlSel(penalty = "SCAD"),
#                        control_outcome = controlOut(penalty = "SCAD"))
# )
#
# expect_equivalent(test_dr_1_scad$output$mean,
#                   0.604803,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_scad$confidence_interval[1] < 0.604803) &
#     (0.604803 < test_dr_1_scad$confidence_interval[2])
# )
#
# # LASSO
# expect_silent(
#   test_dr_1_lasso <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_selection = controlSel(penalty = "lasso"),
#                        control_outcome = controlOut(penalty = "lasso"))
# )
#
# expect_equivalent(test_dr_1_lasso$output$mean,
#                   0.604803,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_lasso$confidence_interval[1] < 0.604803) &
#     (0.604803 < test_dr_1_lasso$confidence_interval[2])
# )
#
# # MCP
# expect_silent(
#   test_dr_1_mcp <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        selection = ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_selection = controlSel(penalty = "MCP"),
#                        control_outcome = controlOut(penalty = "MCP"))
# )
#
# expect_equivalent(test_dr_1_mcp$output$mean,
#                   0.604803,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_dr_1_mcp$confidence_interval[1] < 0.604803) &
#     (0.604803< test_dr_1_mcp$confidence_interval[2])
# )

# MI
# GLM
# test_mi_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "glm")
#
# expect_silent(
#   test_mi_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "glm")
# )
#
# expect_equivalent(test_mi_1$output$mean,
#                   0.6107926,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1$confidence_interval[1] < 0.6107926) &
#     (0.6107926 < test_mi_1$confidence_interval[2])
# )
#
# #NN
# test_mi_1_nn <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "nn")
#
# expect_silent(
#   test_mi_1_nn <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "nn")
# )
#
# expect_equivalent(test_mi_1_nn$output$mean,
#                   0.6410426,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_nn$confidence_interval[1] < 0.6410426) &
#     (0.6410426 < test_mi_1_nn$confidence_interval[2])
# )
#
# # PMM
# test_mi_1_pmm <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "pmm")
#
# expect_silent(
#   test_mi_1_pmm <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "pmm")
# )
#
# expect_equivalent(test_mi_1_pmm$output$mean,
#                   0.6490007,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_pmm$confidence_interval[1] < 0.6490007) &
#     (0.6490007 < test_mi_1_pmm$confidence_interval[2])
# )
#
# test_mi_1_pmm_2 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "pmm",
#                      control_outcome = controlOut(predictive_match = 2))
#
# expect_silent(
#   test_mi_1_pmm_2 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "pmm",
#                        control_outcome = controlOut(predictive_match = 2))
# )
#
# expect_equivalent(test_mi_1_pmm_2$output$mean,
#                   0.7553603,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_pmm_2$confidence_interval[1] < 0.7553603) &
#     (0.7553603 < test_mi_1_pmm_2$confidence_interval[2])
# )

# # SCAD
# test_mi_1_scad <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "glm",
#                      control_inference = controlInf(vars_selection = TRUE),
#                      control_outcome = controlOut(penalty = "SCAD"))
#
# expect_silent(
#   test_mi_1_scad <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "glm",
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_outcome = controlOut(penalty = "SCAD"))
# )
#
# expect_equivalent(test_mi_1_scad$output$mean,
#                   0.6107926,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_scad$confidence_interval[1] < 0.6107926) &
#     (0.6107926 < test_mi_1_scad$confidence_interval[2])
# )
#
# # LASSO
# test_mi_1_lasso <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "glm",
#                      control_inference = controlInf(vars_selection = TRUE),
#                      control_outcome = controlOut(penalty = "lasso"))
#
# expect_silent(
#   test_mi_1_lasso <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "glm",
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_outcome = controlOut(penalty = "lasso"))
# )
#
# expect_equivalent(test_mi_1_lasso$output$mean,
#                   0.6107926,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_lasso$confidence_interval[1] < 0.6107926) &
#     (0.6107926 < test_mi_1_lasso$confidence_interval[2])
# )
#
# # MCP
# test_mi_1_mcp <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                      data = cbop_df_long,
#                      svydesign = popyt_svy,
#                      method_outcome = "glm",
#                      control_inference = controlInf(vars_selection = TRUE),
#                      control_outcome = controlOut(penalty = "MCP"))
#
# expect_silent(
#   test_mi_1_mcp <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
#                        data = cbop_df_long,
#                        svydesign = popyt_svy,
#                        method_outcome = "glm",
#                        control_inference = controlInf(vars_selection = TRUE),
#                        control_outcome = controlOut(penalty = "MCP"))
# )
#
# expect_equivalent(test_mi_1_mcp$output$mean,
#                   0.6107926,
#                   tolerance = 0.1)
#
# expect_true(
#   (test_mi_1_mcp$confidence_interval[1] < 0.6107926) &
#     (0.6107926 < test_mi_1_mcp$confidence_interval[2])
# )
