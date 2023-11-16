# S3methods tests
# test simulate
# set.seed(123)
source_nonprob_p <- read.csv("test1_nonprob.csv")
sample_a <- read.csv("test1_prob.csv")
svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)

test1a <- nonprob(selection = ~ x,
                  target = ~ y1,
                  data = source_nonprob_p,
                  method_selection = "logit",
                  svydesign = svy_a)

expect_silent(
  summary(test1a)
)

expect_silent(
  nobs(test1a)
)

expect_silent(
  pop.size(test1a)
)

expect_silent(
  residuals(test1a)
)

expect_error(
  cooks.distance(test1a)
)

expect_silent(
  hatvalues(test1a)
)

expect_silent(
  logLik(test1a)
)

expect_silent(
  AIC(test1a)
)

expect_silent(
  BIC(test1a)
)

expect_silent(
  confint(test1a)
)

expect_silent(
  vcov(test1a)
)

expect_silent(
  deviance(test1a)
)
