outcomes <- nonprobsvy:::make_outcomes(y1 + y2 ~ x1 + x2)

expect_equal(names(outcomes), c("f", "outcome", "l"))
expect_equal(outcomes$f, c("y1", "y2"))
expect_equal(outcomes$l, 2L)
expect_equal(length(outcomes$outcome), 2L)
expect_equal(stats::formula(outcomes$outcome[[1]]), y1 ~ x1 + x2)
expect_equal(stats::formula(outcomes$outcome[[2]]), y2 ~ x1 + x2)

old_options <- options(warnPartialMatchDollar = TRUE)
on.exit(options(old_options), add = TRUE)

expect_silent(outcomes$outcome[[1]])
