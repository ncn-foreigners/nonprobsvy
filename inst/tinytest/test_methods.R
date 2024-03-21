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

expect_silent(
  residuals(test1a, type = "pearson")
)

expect_silent(
  residuals(test1a, type = "working")
)

expect_silent(
  residuals(test1a, type = "deviance")
)

expect_silent(
  residuals(test1a, "pearsonSTD")
)

# expect_silent(
#   cooks.distance(test1a)
# )
#
# expect_silent(
#   hatvalues(test1a)
# )

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

test2a <- nonprob(outcome = y1 ~ x,
                  data = source_nonprob_p,
                  method_selection = "logit",
                  svydesign = svy_a)
expect_silent(
  summary(test2a)
)

expect_silent(
  nobs(test2a)
)

expect_silent(
  pop.size(test2a)
)

expect_silent(
  residuals(test2a)
)

expect_silent(
  residuals(test2a, type = "pearson")
)

expect_silent(
  residuals(test2a, type = "working")
)

expect_silent(
  residuals(test2a, type = "deviance")
)

expect_silent(
  residuals(test2a, "pearsonSTD")
)

if (isTRUE(tolower(Sys.getenv("TEST_NONPROBSVY_MULTICORE_DEVELOPER")) == "true")) {
  expect_silent(
    cooks.distance(test2a)
  )
  expect_silent(
    hatvalues(test2a)
  )
}

expect_silent(
  logLik(test2a)
)

expect_silent(
  AIC(test2a)
)

expect_silent(
  BIC(test2a)
)

expect_silent(
  confint(test2a)
)

expect_silent(
  vcov(test2a)
)

expect_silent(
  deviance(test2a)
)

test3a <- nonprob(outcome = y1 ~ x,
                  selection = ~ x,
                  data = source_nonprob_p,
                  method_selection = "logit",
                  svydesign = svy_a)
expect_silent(
  summary(test3a)
)

expect_silent(
  nobs(test3a)
)

expect_silent(
  pop.size(test3a)
)

expect_silent(
  residuals(test3a)
)

expect_silent(
  residuals(test3a, type = "pearson")
)

expect_silent(
  residuals(test3a, type = "working")
)

expect_silent(
  residuals(test3a, type = "deviance")
)

expect_silent(
  residuals(test3a, "pearsonSTD")
)

if (isTRUE(tolower(Sys.getenv("TEST_NONPROBSVY_MULTICORE_DEVELOPER")) == "true")) {
  expect_silent(
    cooks.distance(test3a)
  )

  expect_silent(
    hatvalues(test3a)
  )
}

expect_silent(
  logLik(test3a)
)

expect_silent(
  AIC(test3a)
)

expect_silent(
  BIC(test3a)
)

expect_silent(
  confint(test3a)
)

expect_silent(
  vcov(test3a)
)

expect_silent(
  deviance(test3a)
)
