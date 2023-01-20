#' Propensity score estimation by probit model
#' @importFrom maxLik maxLik
#' @importFrom stats pnorm
#' @importFrom stats dnorm
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
      log_like2 <- sum(d_rand * log(1 - invLink2))
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

      t(t(X_nons) %*% (dlink1 / (invLink1 * (1 - invLink1))) - t(X_rand) %*% (d_rand * (dlink2 / (1 - invLink2)))) #jest ok

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

    hess1 <- t(X_nons * ((eta1 * dlink1)/(invLink1 * (1 - invLink1)) - dlink1^2/((invLink1^2) * ((1 - invLink1)^2)) + dlink1/(1 - invLink1)^2)) %*% as.matrix(X_nons)
    hess2 <- - t(X_rand * d_rand * ((eta2 * dlink2)/(1 - invLink2) + (dlink2^2)/((1 - invLink2)^2))) %*% as.matrix(X_rand)
    # jest ok
    hess1 + hess2

    }

  }

  ps_est <- function(X, log_like, gradient, hessian){

    maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                method = "NR", start = rep(0, ncol(X)))

    estim_probit <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    estim_ps <- inv_link(estim_probit %*% t(as.matrix(X)))
    estim_psd <- dlink(estim_probit %*% t(as.matrix(X)))


    return(list("ps" = estim_ps,
                "psd" = estim_psd,
                "grad" = grad,
                "hess" = hess))

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
      linkFun = link,
      linkInv = inv_link,
      linkDer = dlink,
      PropenScore = ps_est
    ),

    class = "method.selection"
  )
}
