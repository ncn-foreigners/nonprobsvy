source("_code_for_all_.R")

set.seed(704)
rcpp_cv_n <- 90
rcpp_cv_p <- 8
rcpp_cv_X <- scale(matrix(rnorm(rcpp_cv_n * rcpp_cv_p), rcpp_cv_n, rcpp_cv_p))
rcpp_cv_eta <- 0.35 * rcpp_cv_X[, 1] - 0.25 * rcpp_cv_X[, 2]
rcpp_cv_R <- rbinom(rcpp_cv_n, size = 1, prob = plogis(rcpp_cv_eta))
rcpp_cv_weights <- rep(1, rcpp_cv_n)
rcpp_cv_pop_totals <- colSums(rcpp_cv_X) + seq_len(rcpp_cv_p) / 100

scad_penalty <- nonprobsvy:::q_lambda_test_cpp(
  par = c(5, 0.5, 1.5, -5),
  lambda = 1,
  penalty = "SCAD",
  a = 3.7
)

expect_equal(scad_penalty[1], 0)
expect_equal(scad_penalty[2], 1)
expect_equal(scad_penalty[3], (3.7 - 1.5) / (3.7 - 1))
expect_equal(scad_penalty[4], 0)

expect_rcpp_cv_shape <- function(result, p, nlambda = 4) {
  expect_equal(length(result$theta_est), p)
  expect_true(all(is.finite(result$theta_est)))
  expect_true(all(result$theta_selected >= 0))
  expect_true(all(result$theta_selected < p))
  expect_true(is.finite(result$lambda))
  expect_equal(length(result$cv_error), nlambda)
  expect_true(all(is.finite(result$cv_error)))
}

expect_rcpp_cv_reproducible <- function(method_selection) {
  set.seed(704)
  cv_first <- nonprobsvy:::cv_nonprobsvy_rcpp(
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

  set.seed(704)
  cv_second <- nonprobsvy:::cv_nonprobsvy_rcpp(
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

  expect_true(identical(cv_first$theta_selected, cv_second$theta_selected))
  expect_true(identical(cv_first$lambda, cv_second$lambda))
  expect_true(identical(cv_first$theta_est, cv_second$theta_est))
  expect_true(identical(cv_first$cv_error, cv_second$cv_error))
}

expect_rcpp_matches_r_gee <- function(method_selection, gee_h_fun) {
  set.seed(42)
  n_rand <- 200
  X_rand <- cbind(1, x1 = rnorm(n_rand), x2 = rnorm(n_rand))
  theta <- c(-0.2, 0.3, -0.25)
  ps <- plogis(drop(X_rand %*% theta))
  X_nons <- X_rand[as.logical(rbinom(n_rand, 1, ps)), , drop = FALSE]
  X_nons <- X_nons[seq_len(min(160, nrow(X_nons))), , drop = FALSE]

  X <- rbind(X_rand, X_nons)
  R <- c(rep(0, nrow(X_rand)), rep(1, nrow(X_nons)))
  weights_rand <- runif(n_rand, 0.8, 1.4)
  weights <- rep(1, nrow(X_nons))
  weights_all <- c(weights_rand, weights)

  r_fit <- suppressWarnings(nonprobsvy:::theta_h_estimation(
    R = R,
    X = X,
    weights_rand = weights_rand,
    weights = weights,
    gee_h_fun = gee_h_fun,
    method_selection = method_selection,
    maxit = 500,
    nleqslv_method = "Broyden",
    nleqslv_global = "dbldog",
    nleqslv_xscalm = "fixed",
    start = rep(0, ncol(X))
  ))

  cpp_fit <- nonprobsvy:::cv_nonprobsvy_rcpp(
    X = unname(X),
    R = R,
    weights_X = weights_all,
    method_selection = method_selection,
    gee_h_fun = gee_h_fun,
    maxit = 500,
    eps = 1e-8,
    lambda_min = 0.01,
    nlambda = 2,
    nfolds = 2,
    penalty = "SCAD",
    a = 3.7,
    pop_totals = NULL,
    verbose = FALSE,
    lambda = 0
  )

  expect_equal(as.vector(cpp_fit$theta_est), as.vector(r_fit$theta_h), tolerance = 1e-5)
}

for (method_selection in c("logit", "probit", "cloglog")) {
  expect_rcpp_cv_reproducible(method_selection)
  for (gee_h_fun in c(1, 2)) {
    expect_rcpp_matches_r_gee(method_selection, gee_h_fun)
  }

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
