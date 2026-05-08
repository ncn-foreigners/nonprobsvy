source("_code_for_all_.R")

set.seed(704)
rcpp_cv_n <- 90
rcpp_cv_p <- 8
rcpp_cv_X <- scale(matrix(rnorm(rcpp_cv_n * rcpp_cv_p), rcpp_cv_n, rcpp_cv_p))
rcpp_cv_eta <- 0.35 * rcpp_cv_X[, 1] - 0.25 * rcpp_cv_X[, 2]
rcpp_cv_R <- rbinom(rcpp_cv_n, size = 1, prob = plogis(rcpp_cv_eta))
rcpp_cv_weights <- rep(1, rcpp_cv_n)
rcpp_cv_pop_totals <- colSums(rcpp_cv_X) + seq_len(rcpp_cv_p) / 100

expect_rcpp_cv_shape <- function(result, p, nlambda = 4) {
  expect_equal(length(result$theta_est), p)
  expect_true(all(is.finite(result$theta_est)))
  expect_true(all(result$theta_selected >= 0))
  expect_true(all(result$theta_selected < p))
  expect_true(is.finite(result$lambda))
  expect_equal(length(result$cv_error), nlambda)
  expect_true(all(is.finite(result$cv_error)))
}

for (method_selection in c("logit", "probit", "cloglog")) {
  set.seed(704)
  cv_result <- nonprobsvy:::cv_nonprobsvy_rcpp(
    X = rcpp_cv_X,
    R = rcpp_cv_R,
    weights_X = rcpp_cv_weights,
    method_selection = method_selection,
    gee_h_fun = 1,
    maxit = 40,
    eps = 1e-4,
    lambda_min = 0.01,
    nlambda = 4,
    nfolds = 2,
    penalty = "lasso",
    a = 3,
    pop_totals = NULL,
    verbose = FALSE
  )

  expect_rcpp_cv_shape(cv_result, rcpp_cv_p)

  cv_result_totals <- nonprobsvy:::cv_nonprobsvy_rcpp(
    X = rcpp_cv_X,
    R = rcpp_cv_R,
    weights_X = rcpp_cv_weights,
    method_selection = method_selection,
    gee_h_fun = 2,
    maxit = 40,
    eps = 1e-4,
    lambda_min = 0.01,
    nlambda = 4,
    nfolds = 2,
    penalty = "lasso",
    a = 3,
    pop_totals = rcpp_cv_pop_totals,
    verbose = FALSE
  )

  expect_rcpp_cv_shape(cv_result_totals, rcpp_cv_p)
}
