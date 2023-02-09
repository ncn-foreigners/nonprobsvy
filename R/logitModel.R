#' Propensity score estimation by logistic regression
#'
#' a method for propensity score estimation using logistic regression basing on dependent variables
#'
#' @importFrom maxLik maxLik
#' @importFrom Matrix Matrix
#' @param ... a
#' @export

logit <- function(...) {

  link <- function(x) {log(x / (1 - x))}
  inv_link <- function(x) {exp(x)/(1 + exp(x))}
  dlink <- function(x) {1 / (x**2 - x)}


  #X_B = X_rand, X_A = X_nons, d = weight
  log_like <- function(X_nons, X_rand, d, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta # linear predictor
      eta2 <- as.matrix(X_rand) %*% theta
      invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)

      log_like1 <- sum(log(invLink1 / (1 - invLink1)))
      log_like2 <- sum(d * log(1 - invLink2))
      log_like1 + log_like2

    }
  }

# d is inverse of probability weights

  gradient <-  function(X_rand, X_nons, d, ...) {

    function(theta) {
      eta2 <- as.matrix(X_rand) %*% theta
      invLink2 <- inv_link(eta2)

      t(t(X_nons) %*% matrix(1, nrow = nrow(X_nons), ncol = 1) - t(X_rand) %*% (d*invLink2))

    }

  }

    hessian <- function(X_rand, X_nons, d, ...) {

      function(theta) {

        eta2 <- as.matrix(X_rand) %*% theta
        invLink2 <- inv_link(eta2)

        - t(as.data.frame(X_rand) * (d * invLink2 * (1 - invLink2))) %*% as.matrix(X_rand)
      }


      }


    ps_est <- function(X, log_like, gradient, hessian, start, optim.method) {

      maxLik_an <- maxLik::maxLik(logLik = log_like, grad = gradient, hess = hessian,
                                  method = optim.method, start = start)
      logit_estim <- maxLik_an$estimate
      grad <- maxLik_an$gradient
      hess <- maxLik_an$hessian

      return(list("ps" = inv_link(logit_estim %*% t(as.matrix(X))),
             "grad" = grad,
             "hess" = hess,
             "theta_hat" = logit_estim))

    }


    variance_covariance1 <- function(X, y, mu, ps, N){ # to fix

      v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * (y - mu)^2))
      v1_ <- 1/N^2 * ((1 - ps)/ps * (y - mu)) %*% X
      v_1 <- t(v1_)

      v_2 <- 0
      for(i in 1:nrow(X)){

        suma <- (1 - ps[i]) * X[i,] %*% t(X[i,])
        v_2 <- v_2 + suma

      }

      v_2 <- 1/N^2 * v_2

      v1_vec <- cbind(v11, v1_)
      v2_mx <- cbind(v_1, v_2)
      V1 <- Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)

      return(V1)
    }

    variance_covariance2 <- function(X, eps, ps, b, n, N){

      s <- eps * as.data.frame(X)
      ci <- n/(n-1) * (1 - ps)
      B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/ps))/sum(ci)
      ei <- (s/ps) - B_hat
      db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)

      D <- (1/N^2) * db_var
      D.var <- b %*% D %*% t(b)

      p <- nrow(D) + 1
      V2 <- Matrix::Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)
      V2[2:p,2:p] <- D

      return(V2)
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
          PropenScore = ps_est,
          VarianceCov1 = variance_covariance1,
          VarianceCov2 = variance_covariance2
        ),
        class = "method.selection"
      )


}
