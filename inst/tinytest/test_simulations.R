source("_code_for_all_.R")


# unit-level data ---------------------------------------------------------

### ipw mle ---------------------------------------------------------------------

ipw_mle <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit")

expect_equal(
  ipw_mle$confidence_interval$lower_bound < kim2019_y_true &
  ipw_mle$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
)


### ipw gee ---------------------------------------------------------------------

ipw_gee <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee"))

expect_equal(
  ipw_gee$confidence_interval$lower_bound < kim2019_y_true &
  ipw_gee$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
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

### mi pmm

mi_pmm <- nonprob(outcome = y1 + y2~x1 + x2,
                 svydesign = kim2019_sample_prob,
                 method_outcome = "pmm",
                 data = kim2019_sample_nonprob,
                 control_inference = control_inf(num_boot = 5))

expect_equal(
  mi_pmm$confidence_interval$lower_bound < kim2019_y_true &
  mi_pmm$confidence_interval$upper_bound > kim2019_y_true,
  c(TRUE, TRUE)
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

ipw_mle <- nonprob(
  selection = ~x1 + x2,
  target = ~y1 + y2,
  pop_total = kim2019_totals,
  data = kim2019_sample_nonprob,
  method_selection = "logit")

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


