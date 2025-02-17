#' @title Probit model for weights adjustment
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{probit_model_nonprobsvy} returns all the methods/objects/functions required to estimate the model, assuming a probit link function.
#' @param ... Additional, optional arguments.
#'
#' @return List with selected methods/objects/functions.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom Matrix Matrix
#' @importFrom survey svyrecvar
#' @keywords internal
#' @export
# must be exported to be visible in c++ script, to consider any other option
probit_model_nonprobsvy <- function(...) {
  link <- function(mu) {
    qnorm(mu)
  } # link
  inv_link <- function(eta) {
    pnorm(eta)
  } # inverse link
  dinv_link <- function(eta) {
    dnorm(eta)
  } # first derivative of inverse link
  dlink <- function(mu) 1 / dnorm(qnorm(mu)) # first derivative of link
  inv_link_rev <- function(eta) {
    -dnorm(eta) / pnorm(eta)^2
  } # first derivative of 1/inv_link
  dinv_link_rev <- function(eta) {
    -dnorm(eta) * (eta + dnorm(eta)) / pnorm(eta)^3
  } # second derivative of 1/inv_link

  log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      # weights_sum <- sum(weights, weights_rand)

      log_like1 <- sum(weights * log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(weights_rand * log(1 - invLink2))
      log_like1 + log_like2
    }
  }

  gradient <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta1 <- as.matrix(X_nons) %*% theta # linear predictor
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      dlink1 <- dinv_link(eta1)
      dlink2 <- dinv_link(eta2)
      # weights_sum <- sum(weights, weights_rand)

      t(t(X_nons) %*% (weights * dlink1 / (invLink1 * (1 - invLink1))) - t(X_rand) %*% (weights_rand * (dlink2 / (1 - invLink2))))
    }
  }

  hessian <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      dlink1 <- dinv_link(eta1)
      dlink2 <- dinv_link(eta2)
      # weights_sum <- sum(weights, weights_rand)

      hess1 <- t(as.data.frame(X_nons) * weights * ((-eta1 * dlink1) / (invLink1 * (1 - invLink1)) - dlink1^2 * (1 - 2 * invLink1) / ((invLink1^2) * ((1 - invLink1)^2)))) %*% as.matrix(X_nons)
      hess2 <- t(as.data.frame(X_rand) * weights_rand * ((-eta2 * dlink2) / (1 - invLink2) + dlink2^2 / ((1 - invLink2)^2))) %*% as.matrix(X_rand)
      hess1 - hess2
    }
  }

  variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, gee_h_fun, weights, pop_totals = NULL) {
    # ensure matrix format and get dimensions
    X <- as.matrix(X)
    n <- nrow(X)
    N <- if (is.null(pop_size)) sum(1 / ps) else pop_size

    # get y values based on pop_size
    y_adj <- if (is.null(pop_size)) {
      weights * (y - mu)
    } else {
      weights * y
    }
    y_sq <- y_adj^2

    # base weights calculation
    w1 <- (1 - ps) / ps^2

    if (est_method == "mle" && is.null(pop_totals)) {
      # MLE specific calculations
      v11 <- sum(w1 * y_sq) / N^2
      v1_ <- (psd * y_adj / ps^2) %*% X / N^2

      # matrix calculations
      v_2 <- matrix(0, ncol = ncol(X), nrow = ncol(X))
      for (i in 1:n) {
        v_2i <- (psd[i] / (ps[i]^2 * (1 - ps[i]))) * (X[i, ] %*% t(X[i, ]))
        v_2 <- v_2 + v_2i
      }
      v_2 <- v_2 / N^2
    } else if (est_method == "gee" && gee_h_fun == 1 || !is.null(pop_totals)) {
      # GEE gee_h_fun=1 or pop_totals case
      v11 <- sum(w1 * y_sq) / N^2
      v1_ <- (w1 * y_adj) %*% X / N^2

      v_2 <- matrix(0, ncol = ncol(X), nrow = ncol(X))
      for (i in 1:n) {
        v_2i <- ((1 - ps[i]) / ps[i]) * (X[i, ] %*% t(X[i, ]))
        v_2 <- v_2 + v_2i
      }
      v_2 <- v_2 / N^2
    } else if (est_method == "gee" && gee_h_fun == 2) {
      # GEE gee_h_fun=2 case
      v11 <- sum(w1 * y_sq) / N^2
      v1_ <- ((1 - ps) / ps * y_adj) %*% X / N^2

      v_2 <- matrix(0, ncol = ncol(X), nrow = ncol(X))
      for (i in 1:n) {
        v_2i <- (1 - ps[i]) * (X[i, ] %*% t(X[i, ]))
        v_2 <- v_2 + v_2i
      }
      v_2 <- v_2 / N^2
    }

    # construct final matrix
    v_1 <- t(v1_)
    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)
  }

  variance_covariance2 <- function(X, svydesign, eps, est_method, gee_h_fun, pop_totals, psd, postStrata = NULL) {
    # get total population size
    N <- sum(1 / svydesign$prob)

    if (!is.null(pop_totals)) {
      # case with population totals
      dim_totals <- length(pop_totals)
      cov <- Matrix::Matrix(0, nrow = dim_totals, ncol = dim_totals, sparse = TRUE)
    } else {
      # adjust probabilities based on method
      if (est_method == "mle") {
        svydesign$prob <- svydesign$prob * ((1 - eps) / psd)
      } else if (est_method == "gee" && gee_h_fun == 2) {
        svydesign$prob <- svydesign$prob * eps
      }

      # ensure X is matrix and properly scaled
      X <- as.matrix(X)
      X_scaled <- sweep(X, 1, svydesign$prob, "/")

      # compute covariance matrix
      cov <- survey::svyrecvar(
        x = X_scaled,
        clusters = svydesign$cluster,
        stratas = svydesign$strata,
        fpcs = svydesign$fpc,
        postStrata = postStrata
      ) / N^2
    }

    # create final sparse matrix
    p <- ncol(cov) + 1
    V2 <- Matrix::Matrix(0, nrow = p, ncol = p, sparse = TRUE)
    V2[2:p, 2:p] <- cov
    V2
  }

  b_vec_ipw <- function(y, mu, ps, psd, eta, X, hess, pop_size, weights, verbose) {
    # get hessian inverse
    hess_inv_neg <- tryCatch(
      solve(-hess),
      error = function(e) {
        if (verbose) message("solve() failed, using ginv() instead.")
        MASS::ginv(-hess)
      }
    )

    # prepare common terms
    X <- as.matrix(X)
    w <- psd / ps^2 * weights
    y_adj <- if (is.null(pop_size)) y - mu else y - mu + 1

    # compute b vector with standard matrix mult
    b <- -(w * y_adj) %*% X %*% hess_inv_neg # TODO opposite sign here (?)

    list(b = b)
  }

  b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights, verbose) {
    # get hessian inverse
    hess_inv <- tryCatch(
      solve(hess),
      error = function(e) {
        if (verbose) message("solve() failed, using ginv() instead.")
        MASS::ginv(hess)
      }
    )

    # prepare matrices and weights
    X <- as.matrix(X)
    w <- psd / ps^2 * weights
    resid <- y - y_pred - h_n

    # compute b vector with standard matrix mult
    -(w * resid) %*% X %*% hess_inv
  }

  t_vec <- function(X, ps, psd, b, y_rand, y_nons, N, weights) { # TODO
    as.vector(psd / (1 - ps)) * X %*% t(as.matrix(b)) + y_rand - 1 / N * sum(weights * y_nons)
  }

  var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
    # weighted residuals
    resid_part <- weights * (y - y_pred - h_n) / ps

    # model adjustment
    model_part <- b %*% t(as.matrix(psd / (ps * (1 - ps)) * as.data.frame(X)))

    # final calculation
    1 / N^2 * sum((1 - ps) * (resid_part - model_part)^2)
  }

    list(
      make_link = link,
      make_dlink = dlink,
      make_log_like = log_like,
      make_gradient = gradient,
      make_hessian = hessian,
      make_link_inv = inv_link,
      make_link_inv_der = dinv_link,
      make_link_inv_rev = inv_link_rev,
      make_link_inv_rev_der = dinv_link_rev,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2,
      b_vec_ipw = b_vec_ipw,
      b_vec_dr = b_vec_dr,
      t_vec = t_vec,
      var_nonprob = var_nonprob
    )
}
