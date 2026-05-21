link_values <- function(link, eta) {
  switch(
    link,
    "logit" = {
      ps <- plogis(eta)
      list(ps = ps, psd = ps * (1 - ps), score_factor = ps)
    },
    "probit" = {
      ps <- pnorm(eta)
      psd <- dnorm(eta)
      list(ps = ps, psd = psd, score_factor = psd / (1 - ps))
    },
    "cloglog" = {
      score_factor <- exp(eta)
      ps <- -expm1(-score_factor)
      list(ps = ps, psd = (1 - ps) * score_factor, score_factor = score_factor)
    }
  )
}

manual_mle_v1 <- function(X, y, mu, ps, psd, pop_size = NULL, weights = rep(1, length(y))) {
  X <- as.matrix(X)
  N <- if (is.null(pop_size)) sum(1 / ps) else pop_size
  y_adj <- if (is.null(pop_size)) weights * (y - mu) else weights * y

  v11 <- sum((1 - ps) / ps^2 * y_adj^2) / N^2
  v1_ <- ((psd / ps^2) * y_adj) %*% X / N^2
  v_2 <- crossprod(X, X * (psd^2 / (ps^2 * (1 - ps)))) / N^2

  rbind(cbind(v11, v1_), cbind(t(v1_), v_2))
}

manual_mle_b <- function(X, y, mu, ps, psd, hess, pop_size = NULL,
                         weights = rep(1, length(y))) {
  y_adj <- if (is.null(pop_size)) y - mu else y
  ((psd / ps^2) * weights * y_adj) %*% as.matrix(X) %*% solve(-hess)
}

manual_mle_v2 <- function(X, svydesign, score_factor) {
  N <- sum(1 / svydesign$prob)
  prob_score <- svydesign$prob / score_factor
  survey::svyrecvar(
    x = sweep(as.matrix(X), 1, prob_score, "/"),
    clusters = svydesign$cluster,
    stratas = svydesign$strata,
    fpcs = svydesign$fpc
  ) / N^2
}

manual_gee_v1 <- function(X, y, mu, ps, pop_size = NULL, gee_h_fun = 1,
                          weights = rep(1, length(y)), pop_totals = NULL) {
  X <- as.matrix(X)
  N <- if (is.null(pop_size)) sum(1 / ps) else pop_size
  y_adj <- if (is.null(pop_size)) weights * (y - mu) else weights * y
  h_inv_ps <- gee_h_fun == 1 || !is.null(pop_totals)

  w11 <- (1 - ps) / ps^2
  w1 <- if (h_inv_ps) w11 else (1 - ps) / ps
  w2 <- if (h_inv_ps) w11 else 1 - ps

  v11 <- sum(w11 * y_adj^2) / N^2
  v1_ <- (w1 * y_adj) %*% X / N^2
  v_2 <- crossprod(X, X * w2) / N^2

  rbind(cbind(v11, v1_), cbind(t(v1_), v_2))
}

manual_gee_v2 <- function(X, svydesign, eps, gee_h_fun = 1) {
  N <- sum(1 / svydesign$prob)
  prob_gee <- if (gee_h_fun == 2) svydesign$prob / eps else svydesign$prob
  survey::svyrecvar(
    x = sweep(as.matrix(X), 1, prob_gee, "/"),
    clusters = svydesign$cluster,
    stratas = svydesign$strata,
    fpcs = svydesign$fpc
  ) / N^2
}

manual_gee_b <- function(X, y, mu, ps, psd, hess, pop_size = NULL,
                         weights = rep(1, length(y))) {
  y_adj <- if (is.null(pop_size)) y - mu else y
  (psd / ps^2 * weights * y_adj) %*% as.matrix(X) %*% solve(-hess)
}

X_nons <- cbind(1, x = c(-0.8, -0.2, 0.4, 1.1))
y_nons <- c(2.5, 1.0, 3.0, 4.5)
case_weights <- c(1.0, 1.4, 0.8, 1.2)
mu <- 2.7
eta_nons <- c(-0.7, -0.1, 0.35, 0.9)
hess <- -matrix(c(7, 1.4, 1.4, 4), nrow = 2)
known_n <- 25

prob_data <- data.frame(
  x = c(-1.0, -0.3, 0.2, 0.75, 1.4),
  base_w = c(4, 5, 3, 6, 4)
)
X_prob <- cbind(1, x = prob_data$x)
eta_prob <- c(-0.9, -0.25, 0.1, 0.65, 1.0)
svydesign <- survey::svydesign(ids = ~1, weights = ~base_w, data = prob_data)

for (link in c("logit", "probit", "cloglog")) {
  nons <- link_values(link, eta_nons)
  prob <- link_values(link, eta_prob)
  method <- method_ps(link)

  expect_equal(
    as.matrix(method$variance_covariance1(
      X = X_nons,
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      pop_size = NULL,
      est_method = "mle",
      gee_h_fun = 1,
      weights = case_weights
    )),
    manual_mle_v1(X_nons, y_nons, mu, nons$ps, nons$psd, weights = case_weights),
    tolerance = 1e-12
  )

  expect_equal(
    as.matrix(method$variance_covariance1(
      X = X_nons,
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      pop_size = known_n,
      est_method = "mle",
      gee_h_fun = 1,
      weights = case_weights
    )),
    manual_mle_v1(
      X_nons,
      y_nons,
      mu,
      nons$ps,
      nons$psd,
      pop_size = known_n,
      weights = case_weights
    ),
    tolerance = 1e-12
  )

  expect_equal(
    method$b_vec_ipw(
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      eta = eta_nons,
      X = X_nons,
      hess = hess,
      pop_size = NULL,
      weights = case_weights,
      verbose = FALSE
    )$b,
    manual_mle_b(X_nons, y_nons, mu, nons$ps, nons$psd, hess, weights = case_weights),
    tolerance = 1e-12
  )

  expect_equal(
    method$b_vec_ipw(
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      eta = eta_nons,
      X = X_nons,
      hess = hess,
      pop_size = known_n,
      weights = case_weights,
      verbose = FALSE
    )$b,
    manual_mle_b(
      X_nons,
      y_nons,
      mu,
      nons$ps,
      nons$psd,
      hess,
      pop_size = known_n,
      weights = case_weights
    ),
    tolerance = 1e-12
  )

  v2 <- method$variance_covariance2(
    X = X_prob,
    svydesign = svydesign,
    eps = prob$ps,
    est_method = "mle",
    gee_h_fun = 1,
    pop_totals = NULL,
    psd = prob$psd
  )

  expect_equal(
    unname(as.matrix(v2[-1, -1])),
    unname(manual_mle_v2(X_prob, svydesign, prob$score_factor)),
    tolerance = 1e-12
  )
}

for (link in c("logit", "probit", "cloglog")) {
  nons <- link_values(link, eta_nons)
  prob <- link_values(link, eta_prob)
  method <- method_ps(link)

  for (gee_h_fun in c(1, 2)) {
    expect_equal(
      as.matrix(method$variance_covariance1(
        X = X_nons,
        y = y_nons,
        mu = mu,
        ps = nons$ps,
        psd = nons$psd,
        pop_size = known_n,
        est_method = "gee",
        gee_h_fun = gee_h_fun,
        weights = case_weights
      )),
      manual_gee_v1(
        X_nons,
        y_nons,
        mu,
        nons$ps,
        pop_size = known_n,
        gee_h_fun = gee_h_fun,
        weights = case_weights
      ),
      tolerance = 1e-12
    )

    v2 <- method$variance_covariance2(
      X = X_prob,
      svydesign = svydesign,
      eps = prob$ps,
      est_method = "gee",
      gee_h_fun = gee_h_fun,
      pop_totals = NULL,
      psd = prob$psd
    )

    expect_equal(
      unname(as.matrix(v2[-1, -1])),
      unname(manual_gee_v2(X_prob, svydesign, prob$ps, gee_h_fun)),
      tolerance = 1e-12
    )
  }

  expect_equal(
    method$b_vec_ipw(
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      eta = eta_nons,
      X = X_nons,
      hess = hess,
      pop_size = NULL,
      weights = case_weights,
      verbose = FALSE,
      est_method = "gee",
      gee_h_fun = 1
    )$b,
    manual_gee_b(X_nons, y_nons, mu, nons$ps, nons$psd, hess, weights = case_weights),
    tolerance = 1e-12
  )

  expect_equal(
    method$b_vec_ipw(
      y = y_nons,
      mu = mu,
      ps = nons$ps,
      psd = nons$psd,
      eta = eta_nons,
      X = X_nons,
      hess = hess,
      pop_size = known_n,
      weights = case_weights,
      verbose = FALSE,
      est_method = "gee",
      gee_h_fun = 2
    )$b,
    manual_gee_b(
      X_nons,
      y_nons,
      mu,
      nons$ps,
      nons$psd,
      hess,
      pop_size = known_n,
      weights = case_weights
    ),
    tolerance = 1e-12
  )
}
