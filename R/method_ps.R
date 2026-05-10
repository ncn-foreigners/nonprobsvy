clamp_prob <- function(p, eps = .Machine$double.eps) {
  p[is.na(p) | p == -Inf] <- eps
  p[p == Inf] <- 1 - eps
  p[p < eps] <- eps
  p[p > 1 - eps] <- 1 - eps
  p
}

stable_softplus <- function(eta) {
  out <- eta
  pos <- eta > 0
  out[pos] <- eta[pos] + log1p(exp(-eta[pos]))
  out[!pos] <- log1p(exp(eta[!pos]))
  out
}

stable_cloglog_rate <- function(eta, eps = .Machine$double.eps) {
  exp(pmin(pmax(eta, log(eps)), log(-log(eps))))
}

mle_ipw_variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, weights) {
  X <- as.matrix(X)
  y <- as.vector(y)
  ps <- as.vector(ps)
  psd <- as.vector(psd)
  weights <- as.vector(weights)
  N <- if (is.null(pop_size)) sum(1 / ps) else pop_size
  y_adj <- if (is.null(pop_size)) weights * (y - mu) else weights * y

  score_weight <- psd / ps^2
  information_weight <- psd^2 / (ps^2 * (1 - ps))

  v11 <- sum((1 - ps) / ps^2 * y_adj^2) / N^2
  v1_ <- (score_weight * y_adj) %*% X / N^2
  v_1 <- t(v1_)
  v_2 <- crossprod(X, X * information_weight) / N^2

  Matrix::Matrix(rbind(cbind(v11, v1_), cbind(v_1, v_2)), sparse = TRUE)
}

mle_ipw_b_vec <- function(y, mu, ps, psd, X, hess, pop_size, weights, verbose) {
  hess_inv_neg <- tryCatch(
    chol2inv(chol(-hess)),
    error = function(e) {
      if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
      MASS::ginv(-hess)
    }
  )

  X <- as.matrix(X)
  y <- as.vector(y)
  ps <- as.vector(ps)
  psd <- as.vector(psd)
  weights <- as.vector(weights)
  y_adj <- if (is.null(pop_size)) y - mu else y
  w <- psd / ps^2 * weights
  b <- (w * y_adj) %*% X %*% hess_inv_neg

  list(b = b)
}

gee_ipw_variance_covariance1 <- function(X, y, mu, ps, pop_size, est_method,
                                         gee_h_fun, weights, pop_totals = NULL) {
  X <- as.matrix(X)
  y <- as.vector(y)
  ps <- as.vector(ps)
  weights <- as.vector(weights)
  N <- if (is.null(pop_size)) sum(1 / ps) else pop_size
  y_adj <- if (is.null(pop_size)) weights * (y - mu) else weights * y

  h_inv_ps <- (est_method == "gee" && gee_h_fun == 1) || !is.null(pop_totals)
  h_x <- est_method == "gee" && gee_h_fun == 2 && is.null(pop_totals)
  if (!h_inv_ps && !h_x) {
    stop("Unsupported GEE h-function for IPW variance.")
  }

  w11 <- (1 - ps) / ps^2
  if (h_inv_ps) {
    w1 <- w11
    w2 <- w11
  } else {
    w1 <- (1 - ps) / ps
    w2 <- 1 - ps
  }

  v11 <- sum(w11 * y_adj^2) / N^2
  v1_ <- (w1 * y_adj) %*% X / N^2
  v_2 <- crossprod(X, X * w2) / N^2

  Matrix::Matrix(rbind(cbind(v11, v1_), cbind(t(v1_), v_2)), sparse = TRUE)
}

legacy_ipw_b_vec <- function(y_adj, w, X, hess, verbose) {
  hess_inv_neg <- tryCatch(
    chol2inv(chol(-hess)),
    error = function(e) {
      if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
      MASS::ginv(-hess)
    }
  )

  X <- as.matrix(X)
  b <- (w * y_adj) %*% X %*% hess_inv_neg

  list(b = b)
}

#' @title Propensity Score Model Functions
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description
#' Function to specify the propensity score (PS) model for the inverse probability weighting estimator.
#' This function provides basic functions logistic regression with a given link function (currently
#' we support `logit`, `probit` and `cloglog`) with additional information about the analytic variance estimator of the mean.
#'
#' @description
#' This is a function returns a `list` of functions that refer to specific estimation methods and variance estimators
#' when whether the IPW alone or the DR estimator is applied. The export of this function is mainly because
#' the functions are used in the variable selection algorithms.
#'
#' Functions starting with `make_log_like`, `make_gradient` and `make_hessian` refer to the maximum likelihood estimation
#' as described in the Chen et al. (2020) paper. These functions take into account different link functions defined through
#' the `link` argument.
#'
#' Functions `make_link`, `make_link_inv`, `make_link_der`, `make_link_inv_der`, `make_link_inv_rev`, and `make_link_inv_rev_der`
#' refer to specific link functions and are used in the estimation process.
#'
#' Functions `variance_covariance1` and `variance_covariance2` refer to the variance estimator of the IPW estimator as
#' defined by Chen et al. (2020).
#'
#' Functions `b_vec_ipw`, `b_vec_dr` and `t_vec` are specific functions defined in the Chen et al. (2020) that are used
#' in the variance estimator of the IPW or the DR.
#'
#' Finally, `var_nonprob` is the non-probability component of the DR estimator as defined by Chen et al. (2020).
#'
#'
#'
#' @param link link for the PS model
#' @param ... Additional, optional arguments.
#'
#' @return A `list` of functions and elements for a specific link function with the following entries:
#'
#' \describe{
#'   \item{make_log_like}{log-likelihood function for a specific link function}
#'   \item{make_gradient}{gradient of the loglik}
#'   \item{make_hessian}{hessian of the loglik}
#'   \item{make_link}{link function}
#'   \item{make_link_inv}{inverse link function}
#'   \item{make_link_der}{first derivative of the link function}
#'   \item{make_link_inv_der}{first derivative of the the inverse link function}
#'   \item{make_link_inv_rev}{defines 1/inv_link}
#'   \item{make_link_inv_rev_der}{first derivative of 1/inv_link}
#'   \item{variance_covariance1}{for the IPW estimator: variance component for the non-probability sample}
#'   \item{variance_covariance2}{for the IPW estimator: variance component for the probability sample}
#'   \item{b_vec_ipw}{for the IPW estimator: the \eqn{b} function as defined in the Chen et al. (2020, sec. 3.2, eq. (9)-(10); sec 4.1)}
#'   \item{b_vec_dr}{for the DR estimator: the \eqn{b} function as defined in the Chen et al. (2020, sec. 3.3., eq. (14); sec 4.1)}
#'   \item{t_vec}{for the DR estimator: the \eqn{b} function as defined in the Chen et al. (2020, sec. 3.3., eq. (14); sec 4.1)}
#'   \item{var_nonprob}{for the DR estimator: non-probability component of the variance for DR estimator}
#'   \item{link}{name of the selected link function for the PS model (character)}
#'   \item{model}{model type (character)}
#' }
#'
#' @examples
#' # Printing information on the model selected
#'
#' method_ps()
#'
#' # extracting specific field
#'
#' method_ps("cloglog")$make_gradient
#'
#' @importFrom Matrix Matrix
#' @importFrom survey svyrecvar
#' @importFrom MASS ginv
#' @importFrom stats dnorm
#' @importFrom stats plogis
#' @importFrom stats pnorm
#' @importFrom stats qlogis
#' @export
method_ps <- function(link = c("logit", "probit", "cloglog"),
                      ...) {

  link_fun <- match.arg(link)

  cloglog <- function(...) {
    link <- function(mu) {
      mu <- clamp_prob(mu)
      log(-log1p(-mu))
    } # link

    inv_link <- function(eta) {
      clamp_prob(-expm1(-stable_cloglog_rate(eta)))
    } # inverse link

    dlink <- function(mu) {
      mu <- clamp_prob(mu)
      1 / ((mu - 1) * log1p(-mu))
    } # first derivative of link

    dinv_link <- function(eta) {
      exp_eta <- stable_cloglog_rate(eta)
      out <- exp_eta * exp(-exp_eta)
      out[!is.finite(out)] <- 0
      out
    } # first derivative of inverse link

    inv_link_rev <- function(eta) {
      exp_eta <- stable_cloglog_rate(eta)
      ps <- inv_link(eta)
      out <- -exp_eta * exp(-exp_eta) / ps^2
      out[!is.finite(out)] <- 0
      out
    } # first derivative of 1/inv_link

    dinv_link_rev <- function(eta) {
      exp_eta <- stable_cloglog_rate(eta)
      exp_neg_exp_eta <- exp(-exp_eta)
      ps <- inv_link(eta)
      out <- exp_eta * exp_neg_exp_eta * (exp_eta - 1) / ps^2 +
        2 * exp_eta^2 * exp_neg_exp_eta^2 / ps^3
      out[!is.finite(out)] <- 0
      out
    } # second derivative of 1/inv_link
    dinv_link_rev2 <- function(eta) {
      exp_eta <- stable_cloglog_rate(eta)
      out <- exp_eta * exp(-exp_eta) * (1 - exp_eta)
      out[!is.finite(out)] <- 0
      out
    }
    log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        # linear predictors
        eta1 <- drop(as.matrix(X_nons) %*% theta)
        eta2 <- drop(as.matrix(X_rand) %*% theta)

        # compute inverse links directly
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)

        # original formula for numerical stability
        log_like1 <- sum(weights * (log(invLink1) - log1p(-invLink1)))
        log_like2 <- sum(weights_rand * log1p(-invLink2))

        log_like1 + log_like2
      }
    }

    gradient <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- as.matrix(X_nons) %*% theta
        eta2 <- as.matrix(X_rand) %*% theta
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)
        # weights_sum <- sum(weights, weights_rand)

        exp_eta1 <- stable_cloglog_rate(eta1)
        exp_eta2 <- stable_cloglog_rate(eta2)

        t(crossprod(X_nons, weights * exp_eta1 / invLink1) - crossprod(X_rand, weights_rand * exp_eta2))

      }
    }

    hessian <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- drop(as.matrix(X_nons) %*% theta)
        eta2 <- drop(as.matrix(X_rand) %*% theta)
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)
        # weights_sum <- sum(weights, weights_rand)

        # Compute the weight vectors first for clarity
        exp_eta1 <- stable_cloglog_rate(eta1)
        exp_eta2 <- stable_cloglog_rate(eta2)

        w_nons <- weights * exp_eta1 / invLink1 * (1 - exp_eta1 / invLink1 + exp_eta1)
        w_rand <- weights_rand * exp_eta2
        mat1 <- crossprod(as.matrix(X_nons), diag(w_nons) %*% as.matrix(X_nons))
        mat2 <- crossprod(as.matrix(X_rand), diag(w_rand) %*% as.matrix(X_rand))

        mat1-mat2
        #mat1 <- crossprod(X_nons * (weights * exp(eta1) / (invLink1) * (1 - exp(eta1) / invLink1 + exp(eta1))),
        #                  as.matrix(X_nons))
        #mat2 <- crossprod(X_rand * weights_rand * exp(eta2), as.matrix(X_rand))
        #mat1-mat2

      }
    }

    variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, gee_h_fun, weights, pop_totals = NULL) {
      if (est_method == "mle" && is.null(pop_totals)) {
        psd <- (1 - ps) * (-log1p(-ps))
        return(mle_ipw_variance_covariance1(
          X = X,
          y = y,
          mu = mu,
          ps = ps,
          psd = psd,
          pop_size = pop_size,
          weights = weights
        ))
      }

      gee_ipw_variance_covariance1(
        X = X,
        y = y,
        mu = mu,
        ps = ps,
        pop_size = pop_size,
        est_method = est_method,
        gee_h_fun = gee_h_fun,
        weights = weights,
        pop_totals = pop_totals
      )
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
          # Probability-sample score factor for cloglog MLE is exp(eta).
          score_factor <- -log1p(-eps)
          svydesign$prob <- svydesign$prob / score_factor
        } else if (est_method == "gee" && gee_h_fun == 2) {
          svydesign$prob <- svydesign$prob / eps
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

    b_vec_ipw <- function(y, mu, ps, X, hess, pop_size, weights, verbose, psd = NULL,
                          eta = NULL, est_method = "mle", gee_h_fun = NULL) {
      if (est_method != "mle") {
        y_adj <- if (is.null(pop_size)) y - mu else y
        return(legacy_ipw_b_vec(
          y_adj = y_adj,
          w = (1 - ps) / ps^2 * exp(eta) * weights,
          X = X,
          hess = hess,
          verbose = verbose
        ))
      }

      score_factor <- if (is.null(eta)) -log1p(-ps) else stable_cloglog_rate(eta)
      psd <- (1 - ps) * score_factor
      mle_ipw_b_vec(
        y = y,
        mu = mu,
        ps = ps,
        psd = psd,
        X = X,
        hess = hess,
        pop_size = pop_size,
        weights = weights,
        verbose = verbose
      )
    }

    b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights, verbose) {
      # get hessian inverse
      hess_inv <- tryCatch(
        chol2inv(chol(hess)),
        error = function(e) {
          if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
          MASS::ginv(hess)
        }
      )

      # prepare common terms
      X <- as.matrix(X)
      w <- ((1 - ps) / ps^2) * weights * exp(eta)
      resid <- y - y_pred - h_n

      # compute b vector with standard matrix mult
      (w * resid) %*% X %*% hess_inv
    }

    t_vec <- function(X, ps, psd, b, y_rand, y_nons, N, weights) {
      as.vector(log(1 - ps)) * tcrossprod(X, as.matrix(b)) + y_rand - 1 / N * sum(weights * y_nons)
    }

    var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
      # prepare matrices
      X <- as.matrix(X)
      b <- as.matrix(b)

      # compute log ratios more stably
      log_ratio <- as.vector(log1p(-ps) - log(ps))

      # compute model adjustment - use standard multiplication
      model_adj <- drop(b %*% t(X * rep(log_ratio, each = ncol(X))))

      # compute weighted residuals
      w_resid <- weights * (y - y_pred - h_n) / ps

      # final variance calculation
      sum((1 - ps) * (w_resid - model_adj)^2) / N^2
    }

    list(
      make_log_like = log_like,
      make_gradient = gradient,
      make_hessian = hessian,
      make_link_fun = link,
      make_link_inv = inv_link,
      make_link_der = dlink,
      make_link_inv_der = dinv_link,
      make_link_inv_rev = inv_link_rev,
      make_link_inv_rev_der = dinv_link_rev,
      make_link_inv_rev_der2 = dinv_link_rev2,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2,
      b_vec_ipw = b_vec_ipw,
      b_vec_dr = b_vec_dr,
      t_vec = t_vec,
      var_nonprob = var_nonprob
    )
  }

  logit <- function(...) {
    link <- function(mu) {
      mu <- clamp_prob(mu)
      qlogis(mu)
    } # link
    inv_link <- function(eta) {
      clamp_prob(plogis(eta))
    } # inverse link
    dlink <- function(mu) {
      mu <- clamp_prob(mu)
      1 / (mu * (1 - mu))
    } # first derivative of link
    dinv_link <- function(eta) {
      p <- inv_link(eta)
      p * (1 - p)
    } # first derivative of inverse link
    inv_link_rev <- function(eta) {
      p <- inv_link(eta)
      -(1 - p) / p
    } # first derivative of 1/inv_link
    dinv_link_rev <- function(eta) {
      p <- inv_link(eta)
      (1 - p) / p
    } # second derivative of 1/inv_link
    dinv_link_rev2 <- function(eta) {
      p <- inv_link(eta)
      p * (1 - p) * (1 - 2 * p)
    }
    log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- drop(X_nons %*% theta) # linear predictor
        eta2 <- drop(X_rand %*% theta)
        log_like1 <- sum(weights * eta1)
        log_like2 <- -sum(weights_rand * stable_softplus(eta2))
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
      if (est_method == "mle" && is.null(pop_totals)) {
        psd <- ps * (1 - ps)
        return(mle_ipw_variance_covariance1(
          X = X,
          y = y,
          mu = mu,
          ps = ps,
          psd = psd,
          pop_size = pop_size,
          weights = weights
        ))
      }

      gee_ipw_variance_covariance1(
        X = X,
        y = y,
        mu = mu,
        ps = ps,
        pop_size = pop_size,
        est_method = est_method,
        gee_h_fun = gee_h_fun,
        weights = weights,
        pop_totals = pop_totals
      )
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

    b_vec_ipw <- function(y, mu, ps, psd, eta, X, hess, pop_size, weights, verbose,
                          est_method = "mle", gee_h_fun = NULL) {
      if (est_method != "mle") {
        y_adj <- if (is.null(pop_size)) y - mu else y
        return(legacy_ipw_b_vec(
          y_adj = y_adj,
          w = (1 - ps) / ps * weights,
          X = X,
          hess = hess,
          verbose = verbose
        ))
      }

      psd <- ps * (1 - ps)
      mle_ipw_b_vec(
        y = y,
        mu = mu,
        ps = ps,
        psd = psd,
        X = X,
        hess = hess,
        pop_size = pop_size,
        weights = weights,
        verbose = verbose
      )
    }

    b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights, verbose) {
      # get hess inverse
      hess_inv <- tryCatch(
        chol2inv(chol(hess)),
        error = function(e) {
          if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
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
      Xb <- tcrossprod(X, as.matrix(b)) # ensure b is matrix
      ps_vec <- as.vector(ps) # ensure ps is vector
      Xb_scaled <- ps_vec * Xb # element-wise mult

      # final calculation
      drop(Xb_scaled + y_rand - mean_nons)
    }

    var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
      # get weighted residuals
      resid <- weights * (y - y_pred - h_n)

      # diff between weighted residuals and model correction
      adj_diff <- resid / ps - drop(tcrossprod(b, X))

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
      make_link_inv_rev_der2 = dinv_link_rev2,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2,
      b_vec_ipw = b_vec_ipw,
      b_vec_dr = b_vec_dr,
      t_vec = t_vec,
      var_nonprob = var_nonprob
    )
  }

  probit <- function(...) {
    link <- function(mu) {
      mu <- clamp_prob(mu)
      qnorm(mu)
    } # link
    inv_link <- function(eta) {
      clamp_prob(pnorm(eta))
    } # inverse link
    dinv_link <- function(eta) {
      dnorm(eta)
    } # first derivative of inverse link
    dlink <- function(mu) 1 / dnorm(qnorm(clamp_prob(mu))) # first derivative of link
    inv_link_rev <- function(eta) {
      -dnorm(eta) / inv_link(eta)^2
    } # first derivative of 1/inv_link
    dinv_link_rev <- function(eta) {
      ps <- inv_link(eta)
      -dnorm(eta) * (eta + dnorm(eta)) / ps^3
    } # second derivative of 1/inv_link
    dinv_link_rev2 <- function(eta) {
      -eta * dnorm(eta)
    }
    log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- as.matrix(X_nons) %*% theta
        eta2 <- as.matrix(X_rand) %*% theta
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)
        # weights_sum <- sum(weights, weights_rand)

        log_like1 <- sum(weights * (log(invLink1) - log1p(-invLink1)))
        log_like2 <- sum(weights_rand * log1p(-invLink2))
        log_like1 + log_like2
      }
    }

    gradient <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- drop(as.matrix(X_nons) %*% theta)
        eta2 <- drop(as.matrix(X_rand) %*% theta)
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)
        dlink1 <- dinv_link(eta1)
        dlink2 <- dinv_link(eta2)
        # weights_sum <- sum(weights, weights_rand)

        w_nons <- weights * dlink1 / (invLink1 * (1 - invLink1))
        w_rand <- weights_rand * (dlink2 / (1 - invLink2))
        t(crossprod(X_nons, w_nons) - crossprod(X_rand, w_rand))

      }
    }

    hessian <- function(X_nons, X_rand, weights, weights_rand, ...) {
      function(theta) {
        eta1 <- drop(as.matrix(X_nons) %*% theta)
        eta2 <- drop(as.matrix(X_rand) %*% theta)
        invLink1 <- inv_link(eta1)
        invLink2 <- inv_link(eta2)
        dlink1 <- dinv_link(eta1)
        dlink2 <- dinv_link(eta2)

        w_nons <- weights * (
          (-eta1 * dlink1) / (invLink1 * (1 - invLink1)) -
            dlink1^2 * (1 - 2 * invLink1) / ((invLink1^2) * ((1 - invLink1)^2))
        )

        w_rand <- weights_rand * (
          (-eta2 * dlink2) / (1 - invLink2) +
            dlink2^2 / ((1 - invLink2)^2)
        )

        crossprod(X_nons, X_nons * w_nons) - crossprod(X_rand,  X_rand * w_rand)

      }
    }

    variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, gee_h_fun, weights, pop_totals = NULL) {
      if (est_method == "mle" && is.null(pop_totals)) {
        return(mle_ipw_variance_covariance1(
          X = X,
          y = y,
          mu = mu,
          ps = ps,
          psd = psd,
          pop_size = pop_size,
          weights = weights
        ))
      }

      gee_ipw_variance_covariance1(
        X = X,
        y = y,
        mu = mu,
        ps = ps,
        pop_size = pop_size,
        est_method = est_method,
        gee_h_fun = gee_h_fun,
        weights = weights,
        pop_totals = pop_totals
      )
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
          svydesign$prob <- svydesign$prob / eps
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

    b_vec_ipw <- function(y, mu, ps, psd, eta, X, hess, pop_size, weights, verbose,
                          est_method = "mle", gee_h_fun = NULL) {
      if (est_method != "mle") {
        y_adj <- if (is.null(pop_size)) y - mu else y
        return(legacy_ipw_b_vec(
          y_adj = y_adj,
          w = psd / ps^2 * weights,
          X = X,
          hess = hess,
          verbose = verbose
        ))
      }

      mle_ipw_b_vec(
        y = y,
        mu = mu,
        ps = ps,
        psd = psd,
        X = X,
        hess = hess,
        pop_size = pop_size,
        weights = weights,
        verbose = verbose
      )
    }

    b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights, verbose) {
      # get hessian inverse
      hess_inv <- tryCatch(
        chol2inv(chol(hess)),
        error = function(e) {
          if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
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
      as.vector(psd / (1 - ps)) * tcrossprod(X, as.matrix(b)) + y_rand - 1 / N * sum(weights * y_nons)
    }

    var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
      # weighted residuals
      resid_part <- weights * (y - y_pred - h_n) / ps

      # model adjustment
      model_part <- tcrossprod(b, as.matrix(psd / (ps * (1 - ps)) * as.data.frame(X)))

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
      make_link_inv_rev_der2 = dinv_link_rev2,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2,
      b_vec_ipw = b_vec_ipw,
      b_vec_dr = b_vec_dr,
      t_vec = t_vec,
      var_nonprob = var_nonprob
    )
  }

  link_funs <- switch(link_fun,
                      "logit" = logit(),
                      "probit" = probit(),
                      "cloglog" = cloglog())

  link_funs$link <- link_fun
  link_funs$model <- "ps"
  class(link_funs) <- "nonprob_method"
  return(link_funs)
}
