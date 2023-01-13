cloglogPs <- function(...){

  link <- function(x) {log(-log(1 - x))}
  inv_link <- function(x) {1 - exp(-exp(x))}
  dlink <- function(x) {1 / (x - 1) * log(1 - x)}

  log_like <- function(X_rand, X_nons, weight = 1, ...){

    function(theta){

      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(weight * log(1 - invLink2))
      log_like1 + log_like2

    }
  }


  gradient <- function(X_rand, X_nons, weight = 1, ...){

    function(theta){

      eta1 <- as.matrix(X_nons)%*%theta
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      grad <- t(t(X_nons) %*% (exp(eta1)/invLink1) - t(X_rand) %*% (d * exp(eta2))) #jest ok
      grad
    }

  }

  hessian <- function(X_rand, X_nons, weight = 1, ...){

    function(theta){

      eta1 <- as.matrix(X_nons)%*%theta
      eta2 <- as.matrix(X_rand) %*%theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      hess <- t(as.data.frame(X_nons) * (exp(eta1)/(invLink1) * (1 - exp(eta1)/invLink1 + exp(eta1)))) %*% as.matrix(X_nons) - t(as.data.frame(X_rand) * weight * exp(eta2)) %*% as.matrix(X_rand) # jest ok
      return(hess)

    }
  }


  ps_est <- function(){

  }

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
