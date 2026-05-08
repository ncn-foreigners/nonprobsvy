local({
  assign(".nonprobsvy_test_nn2_calls", 0L, envir = .GlobalEnv)
  trace(
    "nn2",
    where = asNamespace("RANN"),
    tracer = quote(assign(
      ".nonprobsvy_test_nn2_calls",
      get(".nonprobsvy_test_nn2_calls", envir = .GlobalEnv) + 1L,
      envir = .GlobalEnv
    )),
    print = FALSE
  )
  on.exit({
    untrace("nn2", where = asNamespace("RANN"))
    rm(".nonprobsvy_test_nn2_calls", envir = .GlobalEnv)
  }, add = TRUE)

  toy_svy <- survey::svydesign(
    ids = ~1,
    weights = ~w,
    data = data.frame(x = c(0.5, 4), w = c(2, 3))
  )

  multi_outcome_nn <- nonprob(
    outcome = y1 + y2 ~ x,
    data = data.frame(x = c(0, 2, 5), y1 = c(1, 4, 9), y2 = c(3, 6, 12)),
    svydesign = toy_svy,
    method_outcome = "nn",
    control_outcome = control_out(k = 2),
    se = FALSE
  )

  expect_equal(get(".nonprobsvy_test_nn2_calls", envir = .GlobalEnv), 2L)
  expect_equal(
    unname(multi_outcome_nn$output$mean),
    c(4.9, 7.2),
    tolerance = 1e-12
  )
})
