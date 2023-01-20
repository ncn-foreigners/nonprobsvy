#' Propensity score estimation by logistic regression
#' @importFrom maxLik maxLik
#' @param ... a
#' @export

logit <- function(...) {

  link <- function(x) {log(x / (1 - x))}
  inv_link <- function(x) {exp(x)/(1 + exp(x))}
  dlink <- function(x) {1 / (x**2 - x)}


  #X_B = X_rand, X_A = X_nons, d = weight
  log_like <- function(X_nons, X_rand, d, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(d_rand * log(1 - invLink2))
      log_like1 + log_like2

    }
  }


  gradient <-  function(X_rand, X_nons, d, ...) {

    function(theta) {
      eta2 <- as.matrix(X_rans) %*% theta
      invLink2 <- inv_link(eta2)

      t(t(X_nons) %*% matrix(1, nrow(X_rand), 1) - t(X_rand) %*% (d_rand*invLink2))

    }

  }

    hessian <- function(X_rand, X_nons, d, ...) {

      function(theta) {

        eta2 <- as.matrix(X_rand) %*% theta
        invLink2 <- inv_link(eta2)

        - t(as.data.frame(X_rand) * (d_rand * invLink2 * (1 - invLink2))) %*% as.matrix(X_rand)#jest ok
      }


      }


    ps_est <- function(X, log_like, gradient, hessian) {

      maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                  method = "NR", start = rep(0, ncol(X)))
      logit_estim <- maxLik_an$estimate
      grad <- maxLik_an$gradient
      hess <- maxLik_an$hessian

      return(list("ps" = inv_link(logit_estim %*% t(as.matrix(X))),
             "grad" = grad,
             "hess" = hess))

    }

    # Move to nonprobIPW, MI, DR functions

   # b_var <- function(){

    #  n_X_nons <- nrow(X_nons)

      # a <- 0
      # for(i in 1:n_X_nons){

        # suma <- (1 - estim_ps_nons[i])/estim_ps_nons[i] * (XY_nons[i,1] - mu_hat) * X_nons[i,]
        # a <- a + suma

      # }
     #a <- ((1 - estim_ps_nons)/estim_ps_nons * (XY_nons[,1] - mu_hat)) %*% as.matrix(X_nons)
    #  b <- solve(hess)
    #  a %*% b


  #  }

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
