#' @importFrom maxLik maxLik
#' @importFrom Matrix Matrix

cloglog <- function(...) {

  link <- function(mu) {log(-log(1 - mu))} # link
  inv_link <- function(eta) {1 - exp(-exp(eta))} # inverse link
  dlink <- function(mu) {1 / ((mu - 1) * log(1 - mu))} # first derivative of link
  dinv_link <- function(eta) {exp(eta - exp(eta))} # first derivative of inverse link
  inv_link_rev <- function(eta) {exp(eta + exp(eta))/(1 - exp(-exp(eta)))^2} # first derivative of 1/inv_link

  log_like <- function(X_nons, X_rand, weights, ...) {

    function(theta) {
      eta1 <- as.matrix(X_nons) %*% theta #linear predictor
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(weights * log(1 - invLink2))
      log_like1 + log_like2
    }
  }


  gradient <- function(X_nons, X_rand, weights, ...) {

    function(theta) {
      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      t(t(X_nons) %*% (exp(eta1)/invLink1) - t(X_rand) %*% (weights * exp(eta2)))
    }
  }


  hessian <-  function(X_nons, X_rand, weights, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      t(as.data.frame(X_nons) * (exp(eta1)/(invLink1) * (1 - exp(eta1)/invLink1 + exp(eta1)))) %*% as.matrix(X_nons) - t(as.data.frame(X_rand) * weights * exp(eta2)) %*% as.matrix(X_rand)

    }
  }


  ps_est <- function(X, log_like, gradient, hessian, start, optim_method) {

    maxLik_an <- maxLik::maxLik(logLik = log_like,
                                grad = gradient,
                                hess = hessian,
                                method = optim_method,
                                start = rep(0, length(start)))

    if (maxLik_an$code %in% c(3:7, 100)) {
      switch (as.character(maxLik_an$code),
              "3" = warning("Sam nie wiem co panu jest"),
              "4" = warning("Maxiteration limit reached in fitting ps_est by maxLik."),
              "5" = stop("Inifinite value of log_like in fitting ps_est by maxLik, error code 5"),
              "6" = stop("Inifinite value of gradient in fitting ps_est by maxLik, error code 6"),
              "7" = stop("Inifinite value of hessian in fitting ps_est by maxLik, error code 7"),
              "100" = stop("Error in fitting ps_est with maxLik, error code 100:: Bad start."),
      )
    }

    theta <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    estim_ps <- inv_link(theta %*% t(as.matrix(X)))


    list(ps = estim_ps,
         grad = grad,
         hess = hess,
         theta_hat = theta)
  }


  variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, h) {

    N <- pop_size
    if (est_method == "mle") {
      if (is.null(N)) {
        N <- sum(1/ps)
        v11 <- 1/N^2 * sum((((1 - ps)/ps^2) * (y - mu)^2))
        v1_ <- - 1/N^2 * ((1 - ps)/ps^2 * log(1 - ps) * (y - mu)) %*% X
        v_1 <- t(v1_)
      } else {
        v11 <- 1/N^2 * sum((((1 - ps)/ps^2) * y^2))
        v1_ <- - 1/N^2 * ((1 - ps)/ps^2 * log(1 - ps) * y) %*% X
        v_1 <- t(v1_)
      }

      v_2 <- 0
      for (i in 1:nrow(X)) {
        v_2i <- (1 - ps[i])/ps[i]^2 * log(1-ps[i])^2 * X[i,] %*% t(X[i,])
        v_2 <- v_2 + v_2i
      }
      v_2 <- 1/N^2 * v_2
    } else if (est_method == "gee" && h == "1") {
      if (is.null(N)) {
        N <- sum(1/ps)
        v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * (y - mu)^2))
        v1_ <- 1/N^2 * ((1 - ps)/ps^2 * (y - mu)) %*% X
        v_1 <- t(v1_)
      } else {
        v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * y^2))
        v1_ <- 1/N^2 * ((1 - ps)/ps * y) %*% X
        v_1 <- t(v1_)
      }

      v_2 <- 0
      for(i in 1:nrow(X)){
        v_2i <- (1 - ps[i])/ps[i] * X[i,] %*% t(X[i,])
        v_2 <- v_2 + v_2i
      }
    } else if (est_method == "gee" && h == "2") {
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
    }

    v_2 <- 1/N^2 * v_2
    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    V1 <- Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)
    V1
  }

  variance_covariance2 <- function(X, eps, ps, psd, n, N, est_method, h) {

    if (est_method == "mle") {
      s <- log(1 - eps) * as.data.frame(X)
      ci <- n/(n-1) * (1 - ps)
      B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
      ei <- (s/ps) - B_hat
      db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)
      #D.var <- b %*% D %*% t(b)
    } else if (est_method == "gee"){
      if (h == "1"){
        s <- as.data.frame(X)
        ci <- n/(n-1) * (1 - ps)
        B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
        ei <- (s/ps) - B_hat
        db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)
      } else if (h == "2") {
        s <- eps * as.data.frame(X)
        ci <- n/(n-1) * (1 - ps)
        B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
        ei <- (s/ps) - B_hat
        db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)
      }
    }

    D <- 1/N^2 * db_var
    p <- nrow(D) + 1
    V2 <- Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)
    V2[2:p,2:p] <- D
    V2
  }

  UTB <- function(X, R, weights, ps, eta_pi, mu_der, res) {

    n <- length(R)
    R_rand <- 1 - R

    #print(summary((1 - ps)/ps^2 * as.vector(exp(eta_pi))))
    #print(summary(ps))

    utb <- c(apply(X * R/ps * mu_der - X * R_rand * weights * mu_der, 2, sum),
             apply(X * R * as.vector(inv_link_rev(eta_pi)) * res, 2, sum))/n ## as.vector(exp(eta))

    utb

  }

  b_vec_ipw <- function(y, mu, ps, psd = NULL, eta = NULL, X, hess, pop_size) {

    hess_inv <- solve(hess)
    if (is.null(pop_size)) {
      b <- - ((1 - ps)/ps^2 * exp(eta) * (y - mu)) %*% X %*% hess_inv
    } else {
      b <- - ((1 - ps)/ps^2 * exp(eta) * y) %*% X %*% hess_inv
    }

    list(b = b,
         hess_inv = hess_inv)
  }

  b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess) {
    hess_inv <- solve(hess)
    (((1 - ps)/ps^2) * (y - y_pred - h_n) * exp(eta)) %*% X %*% hess_inv
  }

  t_vec <- function(X, ps, psd, b, y_rand, y_nons, N) {
    as.vector(log(1 - ps)) * X %*% t(as.matrix(b)) + y_rand - 1/N * sum(y_nons)
  }

  var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N) {
    1/N^2 * sum((1 - ps) * (((y - y_pred - h_n)/ps) - b %*% t(as.matrix(log((1 - ps)/ps) * as.data.frame(X))))^2)
  }

  structure(
    list(
      make_log_like = log_like,
      make_gradient = gradient,
      make_hessian = hessian,
      make_link_fun = link,
      make_link_inv = inv_link,
      make_link_der = dlink,
      make_link_inv_der = dinv_link,
      make_link_inv_rev = inv_link_rev,
      make_propen_score = ps_est,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2,
      UTB = UTB,
      b_vec_ipw = b_vec_ipw,
      b_vec_dr = b_vec_dr,
      t_vec = t_vec,
      var_nonprob = var_nonprob
    ),

    class = "method_selection"
  )

}
