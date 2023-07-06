#' @title - logit model
#' @author - Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description -
#'
#' @param ... -
#'
#' @importFrom maxLik maxLik
#' @importFrom Matrix Matrix
#' @importFrom survey svyrecvar
#' @export

logit <- function(...) {

  link <- function(mu) {log(mu / (1 - mu))} # link
  inv_link <- function(eta) {exp(eta)/(1 + exp(eta))} # inverse link
  dlink <- function(mu) {1 / (mu**2 - mu)} # first derivative of link
  dinv_link <- function(eta) {exp(eta) / ((1 + exp(eta))^2)} # first derivative of inverse link
  inv_link_rev <- function(eta) {-exp(-eta)} # first derivative of 1/inv_link
  dinv_link_rev <- function(eta) {exp(-eta)} # second derivative of 1/inv_link


  log_like <- function(X_nons, X_rand, weights, weights_rand, ...) {

    function(theta) {

      eta1 <- as.matrix(X_nons) %*% theta # linear predictor
      eta2 <- as.matrix(X_rand) %*% theta

      #invLink1 <- inv_link(eta1)
      invLink2 <- inv_link(eta2)
      #weights_sum <- sum(weights, weights_rand)

      log_like1 <- sum(weights * eta1)
      log_like2 <- sum(weights_rand * log(1 - invLink2))
      # Reduced form
      # log_like1 <- sum(weights * eta1)
      # log_like2 <- - sum(weights_rand * log(1 + exp(eta2)))
      log_like1 + log_like2
    }
  }

  gradient <-  function(X_nons, X_rand, weights, weights_rand, ...) {

    function(theta) {
      eta2 <- as.matrix(X_rand) %*% theta
      invLink2 <- inv_link(eta2)
      #weights_sum <- sum(weights, weights_rand)
      t(t(X_nons) %*% as.matrix(weights) - t(X_rand) %*% (weights_rand * invLink2)) # matrix(1, nrow = nrow(X_nons), ncol = 1)
    }
  }

    hessian <- function(X_nons, X_rand, weights, weights_rand, ...) {

      function(theta) {
        eta2 <- as.matrix(X_rand) %*% theta
        invLink2 <- inv_link(eta2)
        #weights_sum <- sum(weights, weights_rand)
       - t(as.data.frame(X_rand) * (weights_rand * invLink2 * (1 - invLink2))) %*% as.matrix(X_rand)
      }
    }

    # TODO error and some other models when svydesign is provided with no weights argument (works with method == "BFGS")
    max_lik <- function(X_nons, X_rand, weights, weights_rand, start, control, ...) {

      log_like <- log_like(X_nons,
                           X_rand,
                           weights,
                           weights_rand)

      gradient <- gradient(X_nons,
                           X_rand,
                           weights,
                           weights_rand)

      hessian <- hessian(X_nons,
                         X_rand,
                         weights,
                         weights_rand)

      # minq <- minqa::bobyqa(par = start, fn = log_like)

      if (control$optimizer == "maxLik") {
        ########### maxLik ##########
        maxLik_an <- maxLik::maxLik(logLik = log_like,
                                    grad = gradient,
                                    hess = hessian,
                                    method = control$maxLik_method,
                                    start = start,
                                    printLevel = control$print_level) # Add printLevel to control

        theta <- maxLik_an$estimate
        grad <- maxLik_an$gradient
        hess <- maxLik_an$hessian
        log_likelihood <- log_like(theta)

         if (maxLik_an$code %in% c(3:7, 100)) {
           switch (as.character(maxLik_an$code),
                   "3" = warning("Warning in fitting selection model with maxLik: probably not converged."),
                   "4" = warning("Max iteration limit reached in fitting selection model by maxLik."),
                   "5" = stop("Inifinite value of log_like in fitting selection model by maxLik, error code 5"),
                   "6" = stop("Inifinite value of gradient in fitting selection model by maxLik, error code 6"),
                   "7" = stop("Inifinite value of hessian in fitting selection model by maxLik, error code 7"),
                   "100" = stop("Error in fitting selection model with maxLik, error code 100:: Bad start."),
           )
         }
      } else if (control$optimizer == "optim") { # TODO add optimParallel for high-dimensional data
      ########### optim ##########
        # start <- rep(0, NCOL(X_nons))

        maxLik_an <- stats::optim(fn = log_like,
                                  gr = gradient,
                                  method = control$optim_method,
                                  par = start,
                                  control = list(fnscale = -1,
                                                 trace = control$trace,
                                                 maxit = control$maxit
                                                 ))
        if (maxLik_an$convergence %in% c(1, 10, 51, 52)) {
          switch (as.character(maxLik_an$convergence),
                  "1" = warning("Warning in fitting selection model with optim: the iteration limit maxit had been reached."),
                  "10" = warning("degeneracy of the Nelder–Mead simplex in fitting selection model by optim."),
                  "51" = warning("warning from the L-BFGS-B when fitting by optim."),
                  "52" = stop("indicates an error from the L-BFGS-B method when fitting by optim.")
          )
        }

        theta <- maxLik_an$par
        log_likelihood <- log_like(theta)
        grad <- gradient(theta)
        hess <- hessian(theta)

      }

      list(log_l = log_likelihood,
           grad = grad,
           hess = hess,
           theta_hat = theta)
      }

    variance_covariance1 <- function(X, y, mu, ps, psd, pop_size, est_method, h, weights) {

      N <- pop_size
      if (est_method == "mle" || (est_method == "gee" && h == "2")) {
        if (is.null(N)) {
          N <- sum(1/ps)
          v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * weights * (y - mu)^2))
          v1_ <- 1/N^2 * ((1 - ps)/ps * weights * (y - mu)) %*% X
          v_1 <- t(v1_)
        } else {
          v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * (weights*y)^2))
          v1_ <- 1/N^2 * ((1 - ps)/ps * weights * y) %*% X
          v_1 <- t(v1_)
        }

        v_2 <- 0
        for(i in 1:nrow(X)){
          v_2i <- (1 - ps[i]) * X[i,] %*% t(X[i,])
          v_2 <- v_2 + v_2i
        }
      } else if (est_method == "gee" && h == "1") {
        if (is.null(N)) {
          N <- sum(1/ps)
          v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * weights * (y - mu)^2)) # TODO
          v1_ <- 1/N^2 * ((1 - ps)/ps^2 * weights * (y - mu)) %*% X
          v_1 <- t(v1_)
        } else {
          v11 <- 1/N^2 * sum(((1 - ps)/ps^2 * (weights*y)^2))
          v1_ <- 1/N^2 * ((1 - ps)/ps * weights * y) %*% X
          v_1 <- t(v1_)
        }

        v_2 <- 0
        for(i in 1:nrow(X)){
          v_2i <- (1 - ps[i]) / ps[i] * X[i,] %*% t(X[i,])
          v_2 <- v_2 + v_2i
        }
      }

      v_2 <- 1/N^2 * v_2
      v1_vec <- cbind(v11, v1_)
      v2_mx <- cbind(v_1, v_2)
      V1 <- Matrix::Matrix(rbind(v1_vec, v2_mx), sparse = TRUE)
      V1
    }

    variance_covariance2 <- function(X, svydesign, eps, est_method, h, pop_totals, psd, postStrata = NULL) { #TODO

      N <- sum(1/svydesign$prob)

      #########################
      if (!is.null(pop_totals)) {
        cov <- Matrix::Matrix(nrow = length(pop_totals), ncol = length(pop_totals), data = 0, sparse = TRUE)
      } else {
        if (est_method == "mle" || (est_method == "gee" && h == "2")) {
          svydesign$prob <- as.vector(1/eps * svydesign$prob)
          #X <- as.matrix(eps * as.data.frame(X))
        }
        if (is.null(postStrata)) {
        cov <- 1/N^2 * survey::svyrecvar(X/svydesign$prob, svydesign$cluster, stratas = svydesign$strata, fpcs = svydesign$fpc)
        } else {
          cov <- 1/N^2 * survey::svyrecvar(X/svydesign$prob, svydesign$cluster, stratas = svydesign$strata, fpcs = svydesign$fpc,
                                   postStrata = postStrata)
        }
      }

      p <- ncol(cov) + 1
      V2 <- Matrix::Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE)
      V2[2:p,2:p] <- cov
      V2
    }


    b_vec_ipw <- function(y, mu, ps, psd, eta, X, hess, pop_size, weights) {

      hess_inv <- solve(hess)
      if (is.null(pop_size)) {
        b <- - ((1 - ps)/ps * weights * (y - mu)) %*% X %*% hess_inv # TODO opposite sign here (?)
      } else {
        b <- - ((1 - ps)/ps * weights * y) %*% X %*% hess_inv # TODO opposite sign here (?)
      }
      list(b = b,
           hess_inv = hess_inv)
    }

    b_vec_dr <- function(ps, psd, eta, y, y_pred, mu, h_n, X, hess, weights) {
      hess_inv <- solve(hess)
      - (((1 - ps)/ps) * weights * (y - y_pred - h_n)) %*% X %*% hess_inv
    }

    t_vec <- function(X, ps, psd, b, y_rand, y_nons, N, weights) {
      as.vector(ps) * X %*% t(as.matrix(b)) + y_rand - 1/N * sum(weights * y_nons)
    }

    var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, weights) {
      1/N^2 * sum((1 - ps) * (weights*(y - y_pred - h_n)/ps - b %*% t(X))^2)
    }

      structure(
        list(
          make_log_like = log_like,
          make_gradient = gradient,
          make_hessian = hessian,
          make_link = link,
          make_link_inv = inv_link,
          make_link_der = dlink,
          make_link_inv_der = dinv_link,
          make_link_inv_rev = inv_link_rev,
          make_link_inv_rev_der = dinv_link_rev,
          make_max_lik = max_lik,
          variance_covariance1 = variance_covariance1,
          variance_covariance2 = variance_covariance2,
          b_vec_ipw = b_vec_ipw,
          b_vec_dr = b_vec_dr,
          t_vec = t_vec,
          var_nonprob = var_nonprob
        ),
        class = "method_selection"
      )

}
