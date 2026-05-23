source("_code_for_all_.R")

ipw_est0 <- nonprob(selection = ~ region + private + nace + size,
                    target = ~ single_shift,
                    svydesign = jvs_svy,
                    data = admin,
                    method_selection = "logit",
                    se = TRUE)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
                    target = ~ single_shift,
                    svydesign = jvs_svy,
                    data = admin,
                    method_selection = "logit",
                    se = FALSE)

set.seed(2025)
ipw_est2 <- nonprob(selection = ~ region + private + nace + size,
                    target = ~ single_shift,
                    svydesign = jvs_svy,
                    data = admin,
                    method_selection = "logit",
                    control_inference=control_inf(var_method = "bootstrap",
                                                  num_boot = 20)
)


# test update --------------------------------------------------------------


expect_equal(
  update(ipw_est1, se = TRUE)$output$SE,
  0.00984765550279513
)

expect_equal(
  update(ipw_est1, method_selection = "probit")$output$mean,
  0.708804223550342
)

expect_equal(
  update(ipw_est2, control_inference=control_inf())$output$SE,
  0.00984765550279513
)

# test confint ------------------------------------------------------------

expect_equal(
  confint(ipw_est0),
  structure(list(target = "single_shift", lower_bound = 0.689021847522528,
                 upper_bound = 0.7276239477578), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est0, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.68295701802481,
                 upper_bound = 0.733688777255518), row.names = 1L, class = "data.frame")
)


expect_equal(
  confint(ipw_est1),
  structure(list(target = "single_shift", lower_bound = 0.689021847522528,
                 upper_bound = 0.7276239477578), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est1, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.68295701802481,
                 upper_bound = 0.733688777255518), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est2, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.689343468874017,
                 upper_bound = 0.726465291453178), row.names = 1L, class = "data.frame")
)

# regression test (M7): confint() must honor the requested `level` regardless of
# the fit's control_inf(alpha). A fit with alpha = 0.1 stores a 90% interval;
# confint(level = 0.95) must return the wider 95% interval, not the stored 90% one.
ipw_est_a10 <- nonprob(selection = ~ region + private + nace + size,
                       target = ~ single_shift,
                       svydesign = jvs_svy,
                       data = admin,
                       method_selection = "logit",
                       control_inference = control_inf(alpha = 0.1))
ci95_a10 <- confint(ipw_est_a10, level = 0.95)
expect_equal(
  ci95_a10$upper_bound - ci95_a10$lower_bound,
  2 * stats::qnorm(0.975) * ipw_est_a10$output$SE
)
expect_true(
  (ci95_a10$upper_bound - ci95_a10$lower_bound) >
    (ipw_est_a10$confidence_interval$upper_bound - ipw_est_a10$confidence_interval$lower_bound)
)

