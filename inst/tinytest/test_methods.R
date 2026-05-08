source("_code_for_all_.R")

ipw_est0 <- nonprob(selection = ~ region + private + nace + size,
                    target = ~ single_shift,
                    svydesign = jvs_svy,
                    data = admin,
                    method_selection = "logit",
                    se = T)

ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
                    target = ~ single_shift,
                    svydesign = jvs_svy,
                    data = admin,
                    method_selection = "logit",
                    se = F)

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
  update(ipw_est1, se = T)$output$SE,
  0.00943690366325884
)

expect_equal(
  update(ipw_est1, method_selection = "probit")$output$mean,
  0.708804223550342
)

expect_equal(
  update(ipw_est2, control_inference=control_inf())$output$SE,
  0.00943690366325884
)

# test confint ------------------------------------------------------------

expect_equal(
  confint(ipw_est0),
  structure(list(target = "single_shift", lower_bound = 0.689826906334602,
                 upper_bound = 0.726818888945725), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est0, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.684015044649574,
                 upper_bound = 0.732630750630754), row.names = 1L, class = "data.frame")
)


expect_equal(
  confint(ipw_est1),
  structure(list(target = "single_shift", lower_bound = 0.689826906334602,
                 upper_bound = 0.726818888945725), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est1, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.684015044649574,
                 upper_bound = 0.732630750630754), row.names = 1L, class = "data.frame")
)

expect_equal(
  confint(ipw_est2, level = 0.99),
  structure(list(target = "single_shift", lower_bound = 0.689343468874017,
                 upper_bound = 0.726465291453178), row.names = 1L, class = "data.frame")
)


