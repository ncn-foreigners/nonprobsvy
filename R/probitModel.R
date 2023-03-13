#' Propensity score estimation by probit model
#'
#' a method for propensity score estimation using probit model basing on dependent variables
#'
#' @importFrom maxLik maxLik
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom Matrix Matrix
#' @param ... a
#' @export

probit <- function(...) {

  inv_link <- function(x) {pnorm(x)}
  dinv_link <- function(x) {dnorm(x)}

  log_like <- function(X_nons, X_rand, weights, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta
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

      eta1 <- as.matrix(X_nons)%*%theta #linear predictor
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      dlink1 <- dinv_link(eta1)
      dlink2 <- dinv_link(eta2)

      t(t(X_nons) %*% (dlink1 / (invLink1 * (1 - invLink1))) - t(X_rand) %*% (weights * (dlink2 / (1 - invLink2))))
    }
  }

  hessian <- function(X_nons, X_rand, weights, ...) {

    function(theta) {
    eta1 <- as.matrix(X_nons) %*% theta
    eta2 <- as.matrix(X_rand) %*% theta
    invLink1 <- inv_link(eta1)
    invLink2 <- inv_link(eta2)
    dlink1 <- dinv_link(eta1)
    dlink2 <- dinv_link(eta2)

    hess1 <- t(as.data.frame(X_nons) * ((eta1 * dlink1)/(invLink1 * (1 - invLink1)) - dlink1^2/((invLink1^2) * ((1 - invLink1)^2)) + 2*dlink1/(invLink1*(1 - invLink1)^2))) %*% as.matrix(X_nons)
    hess2 <- t(as.data.frame(X_rand) * weights * ((eta2 * dlink2)/(1 - invLink2) + dlink2^2/((1 - invLink2)^2))) %*% as.matrix(X_rand)

    hess1 - hess2

    }
  }

  ps_est <- function(X, log_like, gradient, hessian, start, optim_method) {

    maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                method = "BFGS", start = start) # NA in gradient for Newton-Raphson method

    estim_probit <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    estim_ps <- inv_link(estim_probit %*% t(X))
    estim_psd <- dinv_link(estim_probit %*% t(X))


    list(ps = estim_ps,
         psd = estim_psd,
         grad = grad,
         hess = hess,
         theta_hat = estim_probit)

  }

  variance_covariance1 <- function(X, y, mu, ps, psd, N = NULL) {

    if (is.null(N)) {

      N <- sum(1/ps)

      v11 <- 1/N^2 * sum((((1 - ps)/ps^2) * (y - mu)^2))
      v1_ <- 1/N^2 *  (psd/ps^2 * (y - mu)) %*% X
      v_1 <- t(v1_)

    } else {

      v11 <- 1/N^2 * sum((((1 - ps)/ps^2) * y^2))
      v1_ <- 1/N^2 *  (psd/ps^2 * y) %*% X
      v_1 <- t(v1_)

    }

    v_2 <- 0
    for (i in 1:nrow(X)) {

      v_2i <- psd[i]/(ps[i]^2 * (1-ps[i])) *  X[i,] %*% t(X[i,])
      v_2 <- v_2 + v_2i

    }

    v_2 <- 1/N^2 * v_2

    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    V1 <- Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)

    V1
  }

  variance_covariance2 <- function(X, eps, ps, psd, n, N) {


    s <- psd/(1-eps) * as.data.frame(X)
    ci <- n/(n-1) * (1 - ps)
    B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
    ei <- (s/ps) - B_hat
    db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)

    D <- (1/N^2) * db_var
    #D.var <- b %*% D %*% t(b) # significantly different than for logit and cloglog

    p <- nrow(D) + 1
    V2 <- Matrix::Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)
    V2[2:p,2:p] <- D

    V2
  }

  # Move to nonprobIPW, MI, DR functions
 # b_var <- function(){

  #  n_A <- nrow(X_nons)

   #  a <- 0
    # for(i in 1:n_A){

     # suma <- estim_psd_nons[i]/estim_ps_nons[i]^2 * (y_nons[i] - mu_hat) * X_nons[i,]
    #  a <- a + suma

    # }

 #   a <- (estim_psd_nons/estim_ps_nons^2 * (y_nons - mu_hat)) %*% as.matrix(X_nons)
  #  b <- solve(hess)
   # a %*% b



  structure(
    list(
      make_log_like = log_like,
      make_gradient = gradient,
      make_hessian = hessian,
      make_link_inv = inv_link,
      make_link_inv_der = dinv_link,
      make_propen_score = ps_est,
      variance_covariance1 = variance_covariance1,
      variance_covariance2 = variance_covariance2
    ),

    class = "method.selection"
  )
}
