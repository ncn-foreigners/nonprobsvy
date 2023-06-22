# TODO:: regenerate data with smarter hiper parameters
# definition of data

# set.seed(123)
# N <- 100000
#
#
# x1 <- rbinom(n = N, size = 3, prob = .4)
# x2 <- rcauchy(n = N)
# x3 <- rgamma(n = N, shape = 2)
# e <- rnorm(n = N)
#
#
# y1 <- 1 + 2 * x1 - .5 * x2 - .4 * x3 + e
# y2 <- rbinom(n = N, size = 1, prob = plogis(-1 + 2 * x2 - x3))
#
# pop <- data.frame(y1, y2, x1, x2, x3, e)
#
# sample_a <- pop[sample(1:N, 500),]
# sample_a$w_a <- N / 500
#
# pp <- exp(-0.5 - 0.1 * x2 - .2 * x3 + .2 * x1) /
#       (1 + exp(-0.5 - 0.1 * x2 - .2 * x3 + .2 * x1))
#
# source_nonprob_p <- pop[rbinom(n = N, size = 1, prob = pp), ]

source_nonprob_p <- read.csv("test1_nonprob.csv")
sample_a <- read.csv("test1_prob.csv")

library(survey)
# IPW ####
expect_silent(
  test1a <- nonprob(
    selection = ~ x1 + x2 + x3,
    target = ~ y1,
    data = source_nonprob_p,
    svydesign = svydesign(
      ids =  ~ 1,
      weights = ~ w_a,
      data = sample_a
    )
  )
)

expect_equivalent(
  test1a$output$mean,
  2.373274,
  tolerance = .01
)

expect_true(
  (test1a$confidence_interval[1] < 2.373274) &
  (2.373274 < test1a$confidence_interval[2])
)

expect_silent(
  test1b <- nonprob(
    selection = ~ x1 + x2 + x3,
    target = ~ y2,
    family_outcome = "binomial",
    data = source_nonprob_p,
    svydesign = svydesign(
      ids =  ~ 1,
      weights = ~ w_a,
      data = sample_a
    )
  )
)

expect_equivalent(
  test1a$output$mean,
  0.24087,
  tolerance = .01
)

expect_true(
  (test1a$confidence_interval[1] < 0.24087) &
  (0.24087 < test1a$confidence_interval[2])
)

# DR  ####
expect_silent(
  test2a <- nonprob(
    selection = ~ x1 + x2 + x3,
    outcome = y1 ~ x1 + x2 + x3,
    data = source_nonprob_p,
    svydesign = svydesign(
      ids =  ~ 1,
      weights = ~ w_a,
      data = sample_a
    )
  )
)

expect_equivalent(
  test2a$output$mean,
  2.373274,
  tolerance = .01
)

expect_true(
  (test2a$confidence_interval[1] < 2.373274) &
  (2.373274 < test2a$confidence_interval[2])
)

expect_silent(
  test2b <- nonprob(
    selection = ~ x1 + x2 + x3,
    outcome = y2 ~ x1 + x2 + x3,
    family_outcome = "binomial",
    data = source_nonprob_p,
    svydesign = svydesign(
      ids =  ~ 1,
      weights = ~ w_a,
      data = sample_a
    )
  )
)

expect_equivalent(
  test2b$output$mean,
  0.24087,
  tolerance = .01
)

expect_true(
  (test2b$confidence_interval[1] < 0.24087) &
    (0.24087 < test2b$confidence_interval[2])
)

# MI  ####
# expect_silent(
#   test3a <- nonprob(
#     outcome = y1 ~ x1 + x2 + x3,
#     data = source_nonprob_p,
#     svydesign = svydesign(
#       ids =  ~ 1,
#       weights = ~ w_a,
#       data = sample_a
#     )
#   )
# )
#
# expect_equivalent(
#   test3a$output$mean,
#   2.373274,
#   tolerance = .01
# )
#
# expect_true(
#   (test3a$confidence_interval[1] < 2.373274) &
#   (2.373274 < test3a$confidence_interval[2])
# )
#
# expect_silent(
#   test3b <- nonprob(
#     outcome = y2 ~ x1 + x2 + x3,
#     family_outcome = "binomial",
#     data = source_nonprob_p,
#     svydesign = svydesign(
#       ids =  ~ 1,
#       weights = ~ w_a,
#       data = sample_a
#     )
#   )
# )
#
# expect_equivalent(
#   test3b$output$mean,
#   0.24087,
#   tolerance = .01
# )
#
# expect_true(
#   (test3b$confidence_interval[1] < 0.24087) &
#     (0.24087 < test3b$confidence_interval[2])
# )

