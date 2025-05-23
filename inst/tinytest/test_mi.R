source("_code_for_all_.R")

# population data only ----------------------------------------------------------------

expect_equal(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    pop_totals = pop_totals,
    method_outcome = "glm",
    data = admin)$output,
  structure(list(mean = 0.703858701257595, SE = 0.00497911011639031),
            row.names = "single_shift", class = "data.frame")
)

## to verification as the SE is very small on comparison to the gaussian
expect_equal(
  nonprob(
    outcome = single_shift~region + private + nace + size,
    pop_totals = pop_totals,
    method_outcome = "glm",
    family_outcome = "binomial",
    data = admin)$output,
  structure(list(mean = 0.757542385494505, SE = 0.00650869641371614),
            row.names = "single_shift", class = "data.frame")
)


# unit-level data ---------------------------------------------------------

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
