# definition of data
# set.seed(123)

# N <- 100000
# n_a <- 500
# n_b <- 1000
# n_b1 <- 0.7*n_b
# n_b2 <- 0.3*n_b
#
# x <- rnorm(N, 2, 1)
# e <- rnorm(N)
# y1 <- 1 + 2*x + e
# y2 <- 3 + x + e
# y3 <- 2.5 + 0.5*x^2 + e
# strata <- x <= 2
#
# pop <- data.frame(x, y1, y2, y3, strata)
#
# # probability sample
# sample_a <- pop[sample(1:N, n_a),]
# sample_a$w_a <- N/n_a
# svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
#
# pop1 <- subset(pop, strata == TRUE)
# pop2 <- subset(pop, strata == FALSE)
#
# # nonprobability sample
# ssource_nonprob_p <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
#                   pop2[sample(1:nrow(pop2), n_b2), ])
# source_nonprob_p$w_b <- N/n_b

library(survey)
source_nonprob_p <- read.csv("test1_nonprob.csv")
sample_a <- read.csv("test1_prob.csv")
svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)

# IPW ####
# logit #
expect_silent(
  test1a <- nonprob(selection = ~ x,
                    target = ~ y1,
                    data = source_nonprob_p,
                    method_selection = "logit",
                    svydesign = svy_a)
  )

expect_equivalent(
  test1a$output$mean,
  4.986351,
  tolerance = .01
)

expect_true(
  (test1a$confidence_interval[1] < 4.986351) &
  (4.986351 < test1a$confidence_interval[2])
)

# cloglog #
expect_silent(
  test1aclog <- nonprob(selection = ~ x,
                    target = ~ y1,
                    data = source_nonprob_p,
                    method_selection = "cloglog",
                    svydesign = svy_a)
)

expect_equivalent(
  test1aclog$output$mean,
  4.98848,
  tolerance = .01
)

expect_true(
  (test1aclog$confidence_interval[1] < 4.98848) &
    (4.98848 < test1aclog$confidence_interval[2])
)

# probit #
expect_silent(
  test1aprob <- nonprob(selection = ~ x,
                        target = ~ y1,
                        data = source_nonprob_p,
                        method_selection = "probit",
                        svydesign = svy_a)
)

expect_equivalent(
  test1aprob$output$mean,
  5.015216,
  tolerance = .01
)

expect_true(
  (test1aprob$confidence_interval[1] < 5.015216) &
    (5.015216 < test1aprob$confidence_interval[2])
)
# DR  ####
# logit #
test2a <- nonprob(selection = ~ x,
                  outcome = y1 ~ x,
                  data = source_nonprob_p,
                  method_selection = "logit",
                  svydesign = svy_a)
expect_silent(
  test2a <- nonprob(selection = ~ x,
                    outcome = y1 ~ x,
                    data = source_nonprob_p,
                    method_selection = "logit",
                    svydesign = svy_a)
)

expect_equivalent(
  test2a$output$mean,
  5.090652,
  tolerance = .01
)

expect_true(
  (test2a$confidence_interval[1] < 5.090652) &
  (5.090652 < test2a$confidence_interval[2])
)

# cloglog #
test2aclog <- nonprob(selection = ~ x,
                  outcome = y1 ~ x,
                  data = source_nonprob_p,
                  method_selection = "cloglog",
                  svydesign = svy_a)
expect_silent(
  test2aclog <- nonprob(selection = ~ x,
                    outcome = y1 ~ x,
                    data = source_nonprob_p,
                    method_selection = "cloglog",
                    svydesign = svy_a)
)

expect_equivalent(
  test2aclog$output$mean,
  5.090655,
  tolerance = .01
)

expect_true(
  (test2aclog$confidence_interval[1] < 5.090655) &
    (5.090655 < test2aclog$confidence_interval[2])
)

# probit #
test2aprob <- nonprob(selection = ~ x,
                  outcome = y1 ~ x,
                  data = source_nonprob_p,
                  method_selection = "probit",
                  svydesign = svy_a)
expect_silent(
  test2aprob <- nonprob(selection = ~ x,
                    outcome = y1 ~ x,
                    data = source_nonprob_p,
                    method_selection = "probit",
                    svydesign = svy_a)
)
expect_equivalent(
  test2aprob$output$mean,
  5.091384,
  tolerance = .01
)

expect_true(
  (test2aprob$confidence_interval[1] < 5.091384) &
    (5.091384< test2aprob$confidence_interval[2])
)
# MI - glm ####
test3a <- nonprob(outcome = y1 ~ x,
                  data = source_nonprob_p,
                  method_selection = "logit",
                  svydesign = svy_a)
expect_silent(
  test3a <- nonprob(outcome = y1 ~ x,
                    data = source_nonprob_p,
                    method_selection = "logit",
                    svydesign = svy_a)
)

expect_equivalent(
  test3a$output$mean,
  5.088456,
  tolerance = .01
)

expect_true(
  (test3a$confidence_interval[1] < 5.088456) &
    (5.088456 < test3a$confidence_interval[2])
)

# MI - nn ####
test3ann <- nonprob(outcome = y1 ~ x,
                  data = source_nonprob_p,
                  svydesign = svy_a,
                  method_outcome = "nn")
expect_silent(
  test3ann <- nonprob(outcome = y1 ~ x,
                    data = source_nonprob_p,
                    svydesign = svy_a,
                    method_outcome = "nn")
)
expect_equivalent(
  test3ann$output$mean,
  5.072862,
  tolerance = .01
)

expect_true(
  (test3ann$confidence_interval[1] < 5.072862) &
    (5.072862 < test3ann$confidence_interval[2])
)
