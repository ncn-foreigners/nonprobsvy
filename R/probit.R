probitPS <- function(...){

  inv_link <- function(x) {pnorm(x)}
  dlink <- function(x) {dnorm(x)}

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
      dlink1 <- dlink(eta1)
      dlink2 <- dlink(eta2)

      grad <- t(t(X_nons) %*% (dlink1 / (invLink1 * (1 - invLink1))) - t(X_rand) %*% (weight * (dlink2 / (1 - invLink2)))) #jest ok
      grad
    }


  }

  hessian <- function(X_rand, X_nons, weight = 1, ...){

    eta1 <- as.matrix(X_nons) %*% theta
    eta2 <- as.matrix(X_rand) %*% theta
    invLink1 <- inv_link(eta1)
    invLink2 <- inv_link(eta2)
    dlink1 <- dlink(eta1)
    dlink2 <- dlink(eta2)

    hess1 <- t(X_nons * ((eta1 * dlink1)/(invLink1 * (1 - invLink1)) - dlink1^2/((invLink1^2) * ((1 - invLink1)^2)) + dlink1/(1 - invLink1)^2)) %*% as.matrix(X_nons)
    hess2 <- - t(X_rand * d * ((eta2 * dlink2)/(1 - invLink2) + (dlink2^2)/((1 - invLink2)^2))) %*% as.matrix(X_rand)
    # jest ok
    hess1 + hess2

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
