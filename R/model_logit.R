#' @title Logit model for weights adjustment
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{logit_model_nonprobsvy} returns all the methods/objects/functions required to estimate the model, assuming a logit link function.
#'
#' @param ... Additional, optional arguments.
#'
#' @return List with selected methods/objects/functions.
#'
#' @seealso
#'
#' [nonprob()] -- for fitting procedure with non-probability samples.
#'
#' @importFrom Matrix Matrix
#' @importFrom survey svyrecvar
#' @importFrom stats plogis
#' @importFrom stats qlogis
#'
#'
#' @keywords internal
#' @export
# must be exported to be visible in c++ script, to consider any other option
logit_model_nonprobsvy <- function(...) {
  link <- function(mu) {
    qlogis(mu)
  } # link
  inv_link <- function(eta) {
    plogis(eta)
  } # inverse link
  dlink <- function(mu) {
    1 / (mu * (1 - mu))
  } # first derivative of link
  dinv_link <- function(eta) {
    p <- plogis(eta)
    p * (1 - p)
  } # first derivative of inverse link
  inv_link_rev <- function(eta) {
    -exp(-eta)
  } # first derivative of 1/inv_link
  dinv_link_rev <- function(eta) {
    exp(-eta)
  } # second derivative of 1/inv_link

  log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta1 <- drop(X_nons %*% theta) # linear predictor
      eta2 <- drop(X_rand %*% theta)
      log_like1 <- sum(weights * eta1)
      log_like2 <- -sum(weights_rand * log1p(exp(eta2)))
      log_like1 + log_like2
    }
  }

  gradient <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta2 <- drop(X_rand %*% theta)
      p2 <- plogis(eta2)
      drop(crossprod(X_nons, weights) - crossprod(X_rand, weights_rand * p2))
    }
  }

  hessian <- function(X_nons, X_rand, weights, weights_rand, ...) {
    function(theta) {
      eta2 <- drop(X_rand %*% theta)
      p2 <- plogis(eta2)
      w <- weights_rand * p2 * (1 - p2)
      -crossprod(X_rand, w * X_rand)
    }
  }

  variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, gee_h_fun, weights, pop_totals = NULL) {
    N <- if (is.null(pop_size)) sum(1 / ps) else pop_size
    n <- ifelse(is.null(dim(X)), length(X), nrow(X))

    # get y values based on N and mu
    if (is.null(pop_size)) {
      y_adj <- weights * (y - mu)
      y_sq <- y_adj^2
    } else {
      y_adj <- weights * y
      y_sq <- y_adj^2
    }

    # get weights based on method
    w1 <- if (est_method == "gee" && gee_h_fun == 1 || !is.null(pop_totals)) {
      (1 - ps) / ps^2
    } else {
      (1 - ps) / ps
    }

    # calc v11 and v1_
    v11 <- sum(w1 * y_sq) / N^2
    v1_ <- (w1 * y_adj) %*% X / N^2 # use standard matrix mult instead of crossprod
    v_1 <- t(v1_)

    # calc v_2
    w2 <- if (est_method == "gee" && gee_h_fun == 1 || !is.null(pop_totals)) {
      (1 - ps) / ps
    } else {
      (1 - ps)
    }

    # calc v_2 with explicit loop instead of lapply
    v_2 <- matrix(0, ncol = ncol(X), nrow = ncol(X))
    for (i in 1:n) {
      v_2 <- v_2 + w2[i] * (X[i, ] %*% t(X[i, ]))
    }
    v_2 <- v_2 / N^2

    # construct final matrix
    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)
  }

  variance_covariance2 <- function(X, svydesign, eps, est_method, gee_h_fun, pop_totals, psd, postStrata = NULL) {
    # get population size
    N <- sum(1 / svydesign$prob)

    if (!is.null(pop_totals)) {
      # case with population totals - init empty covariance matrix
      dim_totals <- length(pop_totals)
      cov <- Matrix::Matrix(0, dim_totals, dim_totals, sparse = TRUE)
    } else {
      # adjust probabilities for MLE/GEE if needed
      if (est_method == "mle" || (est_method == "gee" && gee_h_fun == 2)) {
        svydesign$prob <- svydesign$prob / eps
      }

      # ensure X is matrix and properly scaled
      X <- as.matrix(X)
      X_scaled <- sweep(X, 1, svydesign$prob, "/")

      # calculate covariance based on postStrata
      cov <- survey::svyrecvar(
        x = X_scaled,
        clusters = svydesign$cluster,
        stratas = svydesign$strata,
        fpcs = svydesign$fpc,
        postStrata = postStrata
      ) / N^2
    }

    # construct final sparse matrix with covariance
    p <- ncol(cov) + 1
    V2 <- Matrix::Matrix(0, p, p, sparse = TRUE)
    V2[2:p, 2:p] <- cov
    V2
  }

  b_vec_ipw <- function(y, mu, ps, psd, eta, X, hess, pop_size, weights, verbose) {
    # try matrix inversion - if fails, use ginv
    hess_inv_neg <- tryCatch(
      solve(-hess),
      error = function(e) {
        if (verbose) message("solve() failed, using ginv() instead.")
        MASS::ginv(-hess)
      }
    )

    # prep the weights and adjusted y
    w <- (1 - ps) / ps * weights
    y_adj <- if (is.null(pop_size)) y - mu else y

    # main calculation (explicit steps for debugging)
    weighted_y <- w * y_adj # element-wise mult
    b <- -(weighted_y %*% X) %*% hess_inv_neg # matrix mult step by step

    list(b = b)
  }

  b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights, verbose) {
    # get hess inverse
    hess_inv <- tryCatch(
      solve(hess),
      error = function(e) {
        if (verbose) message("solve() failed, using ginv() instead.")
        MASS::ginv(hess)
      }
    )

    # calc weights and residuals
    w <- ((1 - ps) / ps) * weights
    resid <- y - y_pred - h_n

    # use standard matrix multiplication
    -(w * resid) %*% X %*% hess_inv
  }

  t_vec <- function(X, ps, psd, b, y_rand, y_nons, N, weights) {
    # calc mean of non-sampled values
    mean_nons <- sum(weights * y_nons) / N

    # matrix mult and scaling
    Xb <- X %*% t(as.matrix(b)) # ensure b is matrix
    ps_vec <- as.vector(ps) # ensure ps is vector
    Xb_scaled <- ps_vec * Xb # element-wise mult

    # final calculation
    drop(Xb_scaled + y_rand - mean_nons)
  }

  var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
    # get weighted residuals
    resid <- weights * (y - y_pred - h_n)

    # diff between weighted residuals and model correction
    adj_diff <- resid / ps - drop(b %*% t(X))

    # final variance calc
    sum((1 - ps) * adj_diff^2) / N^2
  }

    list(
      make_log_like = log_like,
      make_gradient = gradient,
      make_hessian = hessian,
      make_link = link,
      make_link_inv = inv_link,
      make_link_der = dlink,
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
