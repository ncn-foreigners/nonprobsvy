#' @importFrom maxLik maxLik
#' @importFrom Matrix Matrix
#' @export

logit <- function(...) {

  link <- function(x) {log(x / (1 - x))}
  inv_link <- function(x) {exp(x)/(1 + exp(x))}
  dlink <- function(x) {1 / (x**2 - x)}


  log_like <- function(X_nons, X_rand, weights, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta # linear predictor
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(weights * log(1 - invLink2))
      log_like1 + log_like2
    }
  }

  gradient <-  function(X_nons, X_rand, weights, ...) {

    function(theta) {
      eta2 <- as.matrix(X_rand) %*% theta
      invLink2 <- inv_link(eta2)
      t(t(X_nons) %*% matrix(1, nrow = nrow(X_nons), ncol = 1) - t(X_rand) %*% (weights*invLink2))
    }
  }

    hessian <- function(X_nons, X_rand, weights, ...) {

      function(theta) {
        eta2 <- as.matrix(X_rand) %*% theta
        invLink2 <- inv_link(eta2)
        - t(as.data.frame(X_rand) * (weights * invLink2 * (1 - invLink2))) %*% as.matrix(X_rand)
      }
    }


    ps_est <- function(X, log_like, gradient, hessian, start, optim_method) {

      maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                  method = optim_method, start = start)
      logit_estim <- maxLik_an$estimate
      grad <- maxLik_an$gradient
      hess <- maxLik_an$hessian

      list(ps = inv_link(logit_estim %*% t(as.matrix(X))),
           grad = grad,
           hess = hess,
           theta_hat = logit_estim)
    }


    variance_covariance1 <- function(X, y, mu, ps, pop_size) { # fixed

      N <- pop_size
      if (is.null(N)) {
        N <- sum(1/ps)
        v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * (y - mu)^2))
        v1_ <- 1/N^2 * ((1 - ps)/ps * (y - mu)) %*% X
        v_1 <- t(v1_)
      } else {
        v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * y^2))
        v1_ <- 1/N^2 * ((1 - ps)/ps * y) %*% X
        v_1 <- t(v1_)
      }

      v_2 <- 0
      for(i in 1:nrow(X)){
        v_2i <- (1 - ps[i]) * X[i,] %*% t(X[i,])
        v_2 <- v_2 + v_2i
      }

      v_2 <- 1/N^2 * v_2
      v1_vec <- cbind(v11, v1_)
      v2_mx <- cbind(v_1, v_2)
      V1 <- Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)
      V1
    }

    variance_covariance2 <- function(X, eps, ps, n, N) {
      s <- eps * as.data.frame(X)
      ci <- n/(n-1) * (1 - ps)
      B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
      ei <- (s/ps) - B_hat
      db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)

      D <- (1/N^2) * db_var
      #D.var <- b %*% D %*% t(b)

      p <- nrow(D) + 1
      V2 <- Matrix::Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)

      ###################### consider using survey package for D estimator
      #svydesign <- stats::update(svydesign,
      #                           eps = as.vector(eps))
      #svydesign_mean <- survey::svymean(~eps, svydesign)

      #var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
      #D <- var_prob

      V2[2:p,2:p] <- D
      V2
    }

      structure(
        list(
          make_log_like = log_like,
          make_gradient = gradient,
          make_hessian = hessian,
          make_link = link,
          make_link_inv = inv_link,
          make_link_der = dlink,
          make_propen_score = ps_est,
          variance_covariance1 = variance_covariance1,
          variance_covariance2 = variance_covariance2
        ),
        class = "method_selection"
      )


}
