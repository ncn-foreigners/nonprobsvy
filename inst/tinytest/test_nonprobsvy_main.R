# Definition of data
library(survey)
library(dplyr)
data <- read.csv("test_data.csv")
data |>
  mutate(zawod_kod2 = factor(zawod_kod2),
         jedna_zmiana=as.numeric(jedna_zmiana))

## probability sample
data |>
  filter(is.na(id_popyt)) |>
  select(id_jednostki, klasa_pr, sek, woj, zawod_kod2, wolne_miejsca_cbop, jedna_zmiana) -> cbop_df

## nonprobability sample
data |>
  filter(!is.na(id_popyt)) |>
  select(id_jednostki, klasa_pr, sek, woj, zawod_kod2, waga, wolne_miejsca) -> popyt_df

popyt_svy <- svydesign(ids=~1, strata = ~klasa_pr+sek+woj, weights = ~waga*wolne_miejsca,
                       data = popyt_df)
cbop_df_long <- cbop_df[rep(1:nrow(cbop_df),cbop_df$wolne_miejsca_cbop), ]

# IPW
expect_silent(
  test_ipw_1 <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
                      target = ~ jedna_zmiana ,
                      data = cbop_df_long,
                      svydesign = popyt_svy,
                      method_selection = "logit")
)

expect_equivalent(test_ipw_1$output$mean,
                  0.621989, # TODO,
                  tolerance = 0.1)

expect_true(
  (test_ipw_1$confidence_interval[1] < 0.621989) & (0.621989 < test_ipw_1$confidence_interval[2])
)

expect_silent(
  test_ipw_2 <- nonprob(selection = ~ klasa_pr + sek + zawod_kod2 + woj,
                      target = ~ jedna_zmiana,
                      data = cbop_df,
                      svydesign = popyt_svy,
                      weights = cbop_df$wolne_miejsca_cbop)
)

expect_equivalent(test_ipw_2$output$mean,
                  0.621989, # TODO,
                  tolerance = 0.1)

expect_true(
  (test_ipw_2$confidence_interval[1] < 0.621989) & (0.621989 < test_ipw_2$confidence_interval[2])
)

# DR
expect_silent(
  test_dr_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
                     selection = ~ klasa_pr + sek + zawod_kod2 + woj,
                     data = cbop_df_long,
                     svydesign = popyt_svy)
)

expect_equivalent(test_dr_1$output$mean,
                  0.6048031,
                  tolerance = 0.1)

expect_true(
  (test_dr_1$confidence_interval[1] < 0.6048031) &
    (0.6048031< test_dr_1$confidence_interval[2])
)

pop_totals <- c(`(Intercept)` = 294294, klasa_prM = 128940, klasa_prS = 71230)
test_dr_2 <- nonprob(selection = ~ klasa_pr,
                    outcome = jedna_zmiana ~ klasa_pr,
                    data = subset(data, !is.na(id_cbop)),
                    pop_totals = pop_totals,
                    control_inference = controlInf(var_method = "analytic")) # TODO warning to connected to algorithm convergence

expect_equivalent(test_dr_2$output$mean,
                  0.6747754,
                  tolerance = 0.1)

expect_true(
  (test_dr_2$confidence_interval[1] < 0.6747754) &
    (0.6747754 < test_dr_2$confidence_interval[2])
)

# MI
test_mi_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
                     data = cbop_df_long,
                     svydesign = popyt_svy,
                     method_outcome = "glm")

expect_silent(
  test_mi_1 <- nonprob(outcome = jedna_zmiana ~ klasa_pr + sek + zawod_kod2 + woj,
                     data = cbop_df_long,
                     svydesign = popyt_svy,
                     method_outcome = "glm")
)

expect_equivalent(test_mi_1$output$mean,
                  0.6107926,
                  tolerance = 0.1)

expect_true(
  (test_mi_1$confidence_interval[1] < 0.6107926) &
    (0.6107926 < test_mi_1$confidence_interval[2])
)
