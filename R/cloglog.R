#' Propensity score estimation by complementary loglog model
#'
#' a method for propensity score estimation using cloglog model basing on dependent variables
#'
#' @importFrom maxLik maxLik
#' @importFrom Matrix Matrix
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
           "hess" = hess,
           "theta_hat" = cloglog_estim))

  }


  variance_covariance1 <- function(X, y, mu, ps, N, hessian){

    v11 <- 1/N^2 * sum((((1 - ps)/ps) * (y - mu)^2))
    v1_ <- - 1/N^2 * ((1 - ps)/ps * log(1 - ps) * (y - mu)) %*% as.matrix(X)
    v_1 <- t(v1_)

    v_2 <- 0
    for(i in 1:nrow(X)){

      suma <- (1 - ps[i])/ps[i]^2 * log(1-ps[i])^2 * t(as.matrix(X[i,])) %*% as.matrix(X[i,])
      v_2 <- v_2 + suma

    }

    v_2 <- 1/N^2 * v_2

    v1_vec <- cbind(v11, v1_)
    v2_mx <- cbind(v_1, v_2)
    V1 <- Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)

    return(V1)
  }

  variance_covariance2 <- function(X, ps, n, N){

    s <- log(1 - ps) * X
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
      linkDer = dlink,
      PropenScore = ps_est,
      VarianceCov1 <- variance_covariance1,
      VarianceCov2 <- variance_covariance2
    ),

    class = "method.selection"
  )

}
