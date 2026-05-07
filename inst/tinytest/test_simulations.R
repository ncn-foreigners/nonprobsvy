source("_code_for_all_.R")

expect_joint_matches_single_outcomes <- function(joint_fit, single_fits) {
  for (outcome_name in names(single_fits)) {
    expect_equal(
      joint_fit$output[outcome_name, , drop = FALSE],
      single_fits[[outcome_name]]$output
    )

    expect_equal(
      joint_fit$confidence_interval[outcome_name, , drop = FALSE],
      single_fits[[outcome_name]]$confidence_interval
    )
  }
}

# unit-level data ---------------------------------------------------------

### ipw mle ---------------------------------------------------------------------

ipw_mle <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit")

ipw_mle_y1 <- nonprob(
  selection = ~x1 + x2,
  target = ~y1,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit")

ipw_mle_y2 <- nonprob(
  selection = ~x1 + x2,
  target = ~y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit")

expect_equal(
  ipw_mle$confidence_interval$lower_bound < kim2019_y_true &
  ipw_mle$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  ipw_mle$output$mean,
  c(3.01644580705596, 1.63078036928783)
)

expect_equal(
  ipw_mle$output$SE,
  c(0.0712074021031954, 0.0524674811724362)
)

expect_joint_matches_single_outcomes(
  ipw_mle,
  list(y1 = ipw_mle_y1, y2 = ipw_mle_y2)
)



### ipw gee ---------------------------------------------------------------------

suppressWarnings(
  ipw_gee <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"))
)

ipw_gee_y1 <- suppressWarnings(nonprob(
  selection = ~x1 + x2,
  target = ~y1,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee")
))

ipw_gee_y2 <- suppressWarnings(nonprob(
  selection = ~x1 + x2,
  target = ~y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee")
))

expect_equal(
  ipw_gee$confidence_interval$lower_bound < kim2019_y_true &
  ipw_gee$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  ipw_gee$output$mean,
  c(2.99726031732032, 1.62129373905383)
)

expect_equal(
  ipw_gee$output$SE,
  c(0.0632650377814916, 0.0498758624332513)
)

expect_joint_matches_single_outcomes(
  ipw_gee,
  list(y1 = ipw_gee_y1, y2 = ipw_gee_y2)
)


ipw_gee_h <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 2,
                                  nleqslv_method = "Newton"))

ipw_gee_h_y1 <- nonprob(
  selection = ~x1 + x2,
  target = ~y1,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 2,
                                  nleqslv_method = "Newton"))

ipw_gee_h_y2 <- nonprob(
  selection = ~x1 + x2,
  target = ~y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 2,
                                  nleqslv_method = "Newton"))

expect_equal(
  ipw_gee_h$confidence_interval$lower_bound < kim2019_y_true &
    ipw_gee_h$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  ipw_gee_h$output$mean,
  c(3.01644580669817, 1.63078036903502)
)

expect_equal(
  ipw_gee_h$output$SE,
  c(0.0712074020763667, 0.0524674811523142)
)

expect_joint_matches_single_outcomes(
  ipw_gee_h,
  list(y1 = ipw_gee_h_y1, y2 = ipw_gee_h_y2)
)

### mi glm

mi_glm <- nonprob(outcome = y1 + y2~x1 + x2,
                  svydesign = kim2019_sample_prob,
                  method_outcome = "glm",
                  data = kim2019_sample_nonprob)

expect_equal(
  mi_glm$confidence_interval$lower_bound < kim2019_y_true &
  mi_glm$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  mi_glm$output$mean,
  c(2.99691596891578, 1.62091690224175)
)

expect_equal(
  mi_glm$output$SE,
  c(0.0623219374193106, 0.0495262154516198)
)


### mi nn

mi_nn <- nonprob(outcome = y1 + y2~x1 + x2,
                  svydesign = kim2019_sample_prob,
                  method_outcome = "nn",
                  data = kim2019_sample_nonprob)


expect_equal(
  mi_nn$confidence_interval$lower_bound < kim2019_y_true &
    mi_nn$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  mi_nn$output$mean,
  c(2.9975501665046, 1.58200852515197)
)

expect_equal(
  mi_nn$output$SE,
  c(0.0655674273365258, 0.0622025880888819)
)

### mi pmm

mi_pmm <- nonprob(outcome = y1 + y2~x1 + x2,
                 svydesign = kim2019_sample_prob,
                 method_outcome = "pmm",
                 data = kim2019_sample_nonprob)

expect_equal(
  mi_pmm$confidence_interval$lower_bound < kim2019_y_true &
  mi_pmm$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

expect_equal(
  mi_pmm$output$mean,
  c(2.99871375079429, 1.63414822522373)
)

expect_equal(
  mi_pmm$output$SE,
  c(0.0650693852678778, 0.0550763525250931)
)



### mi npar

# mi_npar <- nonprob(outcome = y1 + y2 ~ x1 + x2,
#                   svydesign = kim2019_sample_prob,
#                   method_outcome = "npar",
#                   data = kim2019_sample_nonprob)
#
# expect_equal(
#   mi_npar$confidence_interval$lower_bound < kim2019_y_true &
#     mi_npar$confidence_interval$upper_bound > kim2019_y_true,
#   c(TRUE, TRUE)
# )


# pop level data ----------------------------------------------------------

### ipw mle (is gee) ---------------------------------------------------------------------

suppressWarnings(
ipw_mle <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  pop_total = kim2019_totals,
  data = kim2019_sample_nonprob,
  method_selection = "logit")
)

expect_equal(
  ipw_mle$confidence_interval$lower_bound < kim2019_y_true &
  ipw_mle$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)

### mi glm

mi_glm <- nonprob(outcome = y1 + y2 ~ x1 + x2,
                  pop_totals = kim2019_totals,
                  method_outcome = "glm",
                  data = kim2019_sample_nonprob)

expect_equal(
  mi_glm$confidence_interval$lower_bound < kim2019_y_true &
  mi_glm$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)
