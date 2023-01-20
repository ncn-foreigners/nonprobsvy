#' Propensity score estimation by complementary loglog model
#' @importFrom maxLik maxLik
#' @param ... a
#' @export
cloglog <- function(...){

  link <- function(x) {log(-log(1 - x))}
  inv_link <- function(x) {1 - exp(-exp(x))}
  dlink <- function(x) {1 / (x - 1) * log(1 - x)}

  log_like <- function(X_nons, X_rand, d, ...){

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


  gradient <- function(X_rand, X_nons, d, ...){

    function(theta){

      eta1 <- as.matrix(X_nons)%*%theta
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      t(t(X_nons) %*% (exp(eta1)/invLink1) - t(X_rand) %*% (d_rand * exp(eta2))) #jest ok

    }

    }


  hessian <-  function(X_rand, X_nons, d, ...){

    function(theta){

      eta1 <- as.matrix(X_nons)%*%theta
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      t(as.data.frame(X_nons) * (exp(eta1)/(invLink1) * (1 - exp(eta1)/invLink1 + exp(eta1)))) %*% as.matrix(X_nons) - t(as.data.frame(X_rand) * d_rand * exp(eta2)) %*% as.matrix(X_rand) # jest ok

    }
  }


  ps_est <- function(X, log_like, gradient, hessian){

    maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                method = "NR", start = rep(0, ncol(X)))

    cloglog_estim <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    estim_ps <- inv_link(cloglog_estim %*% t(as.matrix(X)))

    return(list("ps" = estim_ps,
           "grad" = grad,
           "hess" = hess))

  }

  # Move to nonprobIPW, MI, DR functions

#  b_var <- function(){

 #   n_X_nons <- nrow(X_nons)

  #  N_est_nons <- sum(1/estim_ps_nons)

   # mu_hat <- (1/N_est_nons) * sum(XY_nons[,1]/estim_psA)

    # n_A <- nrow(XA)

    # a <- 0
    # for(i in 1:n_A){

    #  suma <- ((1 - estim_ps_nons[i])/estim_ps_nons[i]^2) * log(1 - estim_ps_nons[i]) * (XY_nons[,1] - mu_hat) * X_nons[i,] # bez kwadadratu
     # a <- a + suma

    # }

   #  a <- ((1 - estim_ps_nons)/estim_ps_nons^2 * log(1 - estim_ps_nons) * (XY_nons[,1] - mu_hat)) %*% as.matrix(X_nons)
  #  b <- solve(hess)
  #  a %*% b


  structure(
    list(
      MakeLogLike = log_like,
      MakeGradient = gradient,
      MakeHessian = hessian,
      linkFun = link,
      linkInv = inv_link,
      linkDer = dlink
    ),

    class = "method.selection"
  )

}
