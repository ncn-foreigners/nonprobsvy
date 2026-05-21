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

## #71: non-gaussian pop_totals analytic variance must not collapse to 0
mi_pois_totals <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  pop_totals = pop_totals,
  method_outcome = "glm",
  family_outcome = "poisson",
  data = admin)
expect_true(is.finite(mi_pois_totals$output$SE) && mi_pois_totals$output$SE > 0)
expect_true(
  mi_pois_totals$confidence_interval$lower_bound <
    mi_pois_totals$confidence_interval$upper_bound
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

# multi-outcome MI fits should agree with one-outcome fits -------------------

mi_fixture_nonprob <- data.frame(
  x1 = seq(-3, 3, length.out = 30),
  x2 = sin(seq(-3, 3, length.out = 30))
)
mi_fixture_nonprob$y1 <- 2 + 0.6 * mi_fixture_nonprob$x1 -
  0.2 * mi_fixture_nonprob$x2 + 0.1 * cos(seq_len(30))
mi_fixture_nonprob$y2 <- -1 + 0.3 * mi_fixture_nonprob$x1 +
  0.5 * mi_fixture_nonprob$x2 + 0.1 * sin(seq_len(30))

mi_fixture_prob <- data.frame(
  x1 = seq(-2.7, 2.7, length.out = 10),
  x2 = sin(seq(-2.7, 2.7, length.out = 10)),
  w = rep(10, 10)
)
mi_fixture_svy <- svydesign(ids = ~1, weights = ~w, data = mi_fixture_prob)

expect_mi_joint_matches_single <- function(method_outcome,
                                           control_outcome = control_out(k = 1)) {
  mi_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    svydesign = mi_fixture_svy,
    data = mi_fixture_nonprob,
    method_outcome = method_outcome,
    control_outcome = control_outcome
  )
  mi_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    svydesign = mi_fixture_svy,
    data = mi_fixture_nonprob,
    method_outcome = method_outcome,
    control_outcome = control_outcome
  )
  mi_joint <- nonprob(
    outcome = y1 + y2 ~ x1 + x2,
    svydesign = mi_fixture_svy,
    data = mi_fixture_nonprob,
    method_outcome = method_outcome,
    control_outcome = control_outcome
  )

  expect_equal(mi_joint$output["y1", , drop = FALSE], mi_y1$output)
  expect_equal(mi_joint$output["y2", , drop = FALSE], mi_y2$output)
  expect_equal(
    mi_joint$confidence_interval["y1", , drop = FALSE],
    mi_y1$confidence_interval
  )
  expect_equal(
    mi_joint$confidence_interval["y2", , drop = FALSE],
    mi_y2$confidence_interval
  )
  expect_equal(mi_joint$SE["y1", , drop = FALSE], mi_y1$SE)
  expect_equal(mi_joint$SE["y2", , drop = FALSE], mi_y2$SE)
  expect_true(all(is.finite(mi_joint$output$SE)))
  expect_true(all(is.finite(as.matrix(mi_joint$confidence_interval))))
  expect_true(all(is.finite(as.matrix(mi_joint$SE))))
}

expect_mi_joint_matches_single("glm")
expect_mi_joint_matches_single("nn")
expect_mi_joint_matches_single(
  "pmm",
  control_outcome = control_out(k = 1, pmm_reg_engine = "glm")
)
expect_mi_joint_matches_single("npar")

# Regression test (M4): the `npar` analytic variance must stay finite and sane
# even when the loess inclusion-propensity proxy returns improper (<= 0 or
# near-zero) fitted values. With this seed one non-probability unit gets a loess
# fitted of about -0.027 (below the 1/sqrt(N) floor), so the clamp in
# method_npar.R activates; without it the 1/pi_hat^2 weight would be spurious or
# unbounded.
set.seed(4)
m4_N <- 1000L
m4_x <- rnorm(m4_N)
m4_R <- rbinom(m4_N, 1, plogis(-1 + 1.5 * m4_x))
m4_y <- 1 + 2 * m4_x + rnorm(m4_N)
m4_nonprob <- data.frame(x = m4_x, y = m4_y)[m4_R == 1, ]
m4_idxB <- sample.int(m4_N, 200)
m4_svy <- svydesign(ids = ~1, weights = ~w,
                    data = data.frame(x = m4_x[m4_idxB], w = m4_N / 200))
m4_npar <- nonprob(outcome = y ~ x, data = m4_nonprob, svydesign = m4_svy,
                   pop_size = m4_N, method_outcome = "npar")
expect_true(is.finite(m4_npar$output$SE) && m4_npar$output$SE > 0 && m4_npar$output$SE < 0.3)
expect_equal(m4_npar$output$SE, 0.147924, tolerance = 1e-4)
