source("_code_for_all_.R")

expect_equal(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin)$output$mean,
  0.7038587
)

expect_equal(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    svydesign = jvs_svy,
    data = admin)$output$mean,
  0.7038587
)
