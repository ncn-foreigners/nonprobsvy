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

probit <- function(...){

  inv_link <- function(x) {pnorm(x)}
  dlink <- function(x) {dnorm(x)}

  log_like <- function(X_rand, X_nons, d, ...){

    function(theta){

      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(d * log(1 - invLink2))
      log_like1 + log_like2

    }
  }

  gradient <-  function(X_rand, X_nons, d, ...){

    function(theta){

      eta1 <- as.matrix(X_nons)%*%theta
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      dlink1 <- dlink(eta1)
      dlink2 <- dlink(eta2)

      t(t(X_nons) %*% (dlink1 / (invLink1 * (1 - invLink1))) - t(X_rand) %*% (d * (dlink2 / (1 - invLink2)))) #jest ok
    }

  }

  hessian <- function(X_rand, X_nons, d, ...){

    function(theta){
    eta1 <- as.matrix(X_nons) %*% theta
    eta2 <- as.matrix(X_rand) %*% theta
    invLink1 <- inv_link(eta1)
    invLink2 <- inv_link(eta2)
    dlink1 <- dlink(eta1)
    dlink2 <- dlink(eta2)

    hess1 <- t(as.data.frame(X_nons) * ((eta1 * dlink1)/(invLink1 * (1 - invLink1)) - dlink1^2/((invLink1^2) * ((1 - invLink1)^2)) + dlink1/(1 - invLink1)^2)) %*% as.matrix(X_nons)
    hess2 <- - t(as.data.frame(X_rand) * d * ((eta2 * dlink2)/(1 - invLink2) + (dlink2^2)/((1 - invLink2)^2))) %*% as.matrix(X_rand)
    # jest ok
    hess1 + hess2

    }

  }

  ps_est <- function(X, log_like, gradient, hessian, start){

    maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                method = "BFGS", start = start)

    estim_probit <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    estim_ps <- inv_link(estim_probit %*% t(X))
    estim_psd <- dlink(estim_probit %*% t(X))


    return(list("ps" = estim_ps,
                "psd" = estim_psd,
                "grad" = grad,
                "hess" = hess,
                "theta_hat" = estim_probit))

  }

  variance_covariance1 <- function(X, y, mu, ps, psd, N){

    v11 <- 1/N^2 * sum((((1 - ps)/ps^2) * (y - mu)^2))
    v1_ <- 1/N^2 *  (psd/ps^2 * (y - mu)) %*% X
    v_1 <- t(v1_)

    v_2 <- 0
    for(i in 1:nrow(X)){

      suma <- psd[i]/(ps[i]^2 * (1-ps[i])) *  X[i,] %*% t(X[i,])
      v_2 <- v_2 + suma

    }

    v_2 <- 1/N^2 * v_2

    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    V1 <- Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)

    return(V1)
  }

  variance_covariance2 <- function(X, eps, ps, psd, b, n, N){

    s <- psd/(1-eps) * as.data.frame(X)
    ci <- n/(n-1) * (1 - ps)
    B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
    ei <- (s/ps) - B_hat
    db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)

    D <- (1/N^2) * db_var
    D.var <- b %*% D %*% t(b)

    p <- nrow(D) + 1
    V2 <- Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)
    V2[2:p,2:p] <- D

    return(V2)
  }

  # Move to nonprobIPW, MI, DR functions
 # b_var <- function(){

  #  n_A <- nrow(X_nons)

   #  a <- 0
    # for(i in 1:n_A){

     # suma <- estim_psd_nons[i]/estim_ps_nons[i]^2 * (y_nons[i] - mu_hat) * X_nons[i,] #bez kwadratu
    #  a <- a + suma

    # }

 #   a <- (estim_psd_nons/estim_ps_nons^2 * (y_nons - mu_hat)) %*% as.matrix(X_nons)
  #  b <- solve(hess)
   # a %*% b



  structure(
    list(
      MakeLogLike = log_like,
      MakeGradient = gradient,
      MakeHessian = hessian,
      linkInv = inv_link,
      linkDer = dlink,
      PropenScore = ps_est,
      VarianceCov1 = variance_covariance1,
      VarianceCov2 = variance_covariance2
    ),

    class = "method.selection"
  )
}
