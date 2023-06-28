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

# DR  ####
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


# MI  ####
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
  test2a$output$mean,
  5.088456,
  tolerance = .01
)

expect_true(
  (test2a$confidence_interval[1] < 5.088456) &
    (5.088456 < test2a$confidence_interval[2])
)
