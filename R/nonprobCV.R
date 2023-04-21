# These functions are only used internally in the nonprobSel function, so there is no need for documenting them
#' @useDynLib nonprobsvy
#' @import RcppArmadillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats glm
#' @importFrom stats residuals
#' @export


cv_nonprobsvy <- function(X,
                          R,
                          weights_X,
                          method_selection,
                          h,
                          maxit,
                          eps,
                          lambda,
                          lambda_min,
                          nlambda,
                          nfolds) {

  if(!is.null(lambda)) {
    lambda <- lambda
  } else {
    loc_nons <- which(R == 1)
    loc_rand <- which(R == 0)
    X_nons <- cbind(X[loc_nons,], weights_X[loc_nons], R[loc_nons])
    X_rand <- cbind(X[loc_rand,], weights_X[loc_rand], R[loc_rand])

    lambdas <- setup_lambda(X = X,
                            y = R,
                            weights = weights_X,
                            method_selection = method_selection,
                            lambda_min = lambda_min,
                            nlambda = nlambda)

    X_nons <- X_nons[sample(nrow(X_nons)), ]
    X_rand <- X_rand[sample(nrow(X_rand)), ]

    folds_nons <- cut(seq(1, nrow(X_nons)), breaks = nfolds, labels = FALSE)
    folds_rand <- cut(seq(1, nrow(X_rand)), breaks = nfolds, labels = FALSE)

    sample_nons <- sample(1:nfolds, nfolds, replace = FALSE)
    sample_rand <- sample(1:nfolds, nfolds, replace = FALSE)

    loss_theta_av <- sapply(lambdas, function(lambda) { # large elapsed time with parallel package
      loss_theta_vec <- sapply(1:nfolds, function(i) {
        # train data for X_nons
        idx_nons <- which(folds_nons == sample_nons[i], arr.ind = TRUE)
        X_nons_train <- X_nons[-idx_nons, ]
        # test data for X_nons
        X_nons_test <- X_nons[idx_nons, ]

        # train data for X_rand
        idx_rand <- which(folds_rand == sample_rand[i], arr.ind = TRUE)
        X_rand_train <- X_rand[-idx_rand, ]
        # test data for X_rand
        X_rand_test <- X_rand[idx_rand, ]

        X_train <- rbind(X_rand_train, X_nons_train)
        X_test <- rbind(X_rand_test, X_nons_test)
        ncols <- ncol(X_test)
        idxx <- 1:(ncols - 2)

        theta_est <- fit_nonprobsvy(X = X_train[, idxx], # fit_nonprobsvy_cpp
                                    R = X_train[, ncols],
                                    weights = X_train[, ncols - 1],
                                    method_selection = method_selection,
                                    h = h,
                                    lambda = lambda,
                                    maxit = maxit,
                                    eps = eps)
            if (any(theta_est == 0)) {
          idxx <- idxx[-which(theta_est == 0)]
        }
        X_testloss <- X_test[,idxx]
        R_testloss <- X_test[, ncols]
        weights_testloss <- X_test[, ncols-1]
        loss_theta(par = theta_est[idxx],
                   X = X_testloss,
                   R = R_testloss,
                   weights = weights_testloss,
                   h = h,
                   method_selection = method_selection)
        print(theta_est)
      })
      mean(loss_theta_vec)}
  )
    min <- which.min(loss_theta_av)
    lambda <- lambdas[min]
    }
  theta_est <- fit_nonprobsvy(X = X,
                              R = R,
                              weights = weights_X,
                              method_selection = method_selection,
                              h = h,
                              lambda = lambda,
                              maxit = maxit,
                              eps = eps,
                              warn = TRUE)
  theta_selected <- which(theta_est != 0) - 1
  list(min = min,
       lambda = lambda,
       theta_est = theta_est,
       theta_selected = theta_selected)
}


# code for the function comes from the ncvreg package
setup_lambda <- function(X,
                         y,
                         weights,
                         method_selection,
                         lambda_min,
                         nlambda,
                         alpha = 1,
                         log_lambda = FALSE,
                         ...) { #consider penalty factor here

  #fit <- glm.fit(x = X, y = y, weights = weights, family = binomial(link = method_selection))
  fit <- stats::glm(y~1,
                    weights = weights,
                    family = binomial(link = method_selection))

  n <- length(y)
  p <- ncol(X)
  w <- fit$weights
  r <- as.matrix(stats::residuals(fit, "working") * w)
  zmax <- max(crossprod(X, r))/n
  lambda_max <- zmax/alpha


  if (log_lambda) { # lambda sequence on log-scale
    if (lambda_min==0) {
      lambda <- c(exp(seq(log(lambda_max), log(.001*lambda_max), length=nlambda-1)), 0)
    } else {
      lambda <- exp(seq(log(lambda_max), log(lambda_min*lambda_max), length=nlambda))
    }
  } else { # lambda sequence on linear-scale
    if (lambda_min==0) {
      lambda <- c(seq(lambda_max, 0.001*lambda_max, length = nlambda-1), 0)
    } else {
      lambda <- seq(lambda_max, lambda_min*lambda_max, length = nlambda)
    }
  }
  lambda
}


fit_nonprobsvy <- function(X,
                           R,
                           weights,
                           method_selection,
                           h,
                           lambda,
                           maxit,
                           eps,
                           warn = FALSE,
                           pop_totals = NULL,
                           ...) {

  p <- ncol(X)
  init_theta <- rep(0, p)
  # variables selection using score equation for theta
  par0 <- init_theta
  LAMBDA <- rep(0, p)
  it <- 0
  for(jj in 1:maxit) {
    it <- it + 1
    if (warn) {
      if (it == maxit) {
        warning("Convergence not obtained in ", maxit, " iterations of fitting algorithm for variables selection", sep="")
        break
      }
    }

    u_theta0 <- u_theta(R = R,
                        X = X,
                        weights = weights,
                        h = h,
                        method_selection = method_selection,
                        pop_totals = pop_totals)

    u_theta0_der <- u_theta_der(R = R,
                                X = X,
                                weights = weights,
                                h = h,
                                method_selection = method_selection,
                                pop_totals = pop_totals)

    LAMBDA <- abs(q_lambda(par0, lambda))/(eps + abs(par0))
    par <- par0 + MASS::ginv(u_theta0_der(par0) + diag(LAMBDA)) %*% (u_theta0(par0) - diag(LAMBDA) %*% par0) # perhaps 'solve' function instead of 'ginv'
    # equation (13) in article
    if (sum(abs(par - par0)) < eps) break;
    if (sum(abs(par - par0)) > 1000) break;

    par0 <- as.vector(par)
  }

  par <- as.vector(par)
  par[abs(par) < 0.001] <- 0
  theta_est <- par
  theta_est
}

# score equation for theta, used in variable selection
u_theta <- function(R,
                    X,
                    weights,
                    method_selection,
                    h,
                    N = NULL,
                    pop_totals = NULL,
                    pop_size = NULL
                    ) {


  method <- get_method(method_selection)
  inv_link <- method$make_link_inv

  function(par) {

    theta <- as.matrix(par)
    n <- length(R)
    X0 <- as.matrix(X)
    eta_pi <- X0 %*% theta
    ps <- inv_link(eta_pi)
    R_rand <- 1 - R
    ps <- as.vector(ps)
    N_nons <- sum(1/ps)

    R <- as.vector(R) # <------ required if Rcpp
    weights <- as.vector(weights) # <------ required if Rcpp
    R_rand <- as.vector(R_rand) # <------ required if Rcpp

    if (is.null(pop_totals)) {
      eq <- switch(h,
                   "1" = c(apply(X0 * R/ps - X0 * R_rand * weights, 2, sum))/N_nons,
                   "2" = c(apply(X0 * R - X0 * R_rand * weights * ps, 2, sum))/N_nons)
    } else {
      eq <- (c(apply(X0 * R/ps, 2, sum)) - c(N_nons, pop_totals))/N_nons
    }
    eq
  }
}


# derivative of score equation for theta, used in variable selection

u_theta_der <-  function(R,
                         X,
                         weights,
                         method_selection,
                         h,
                         N = NULL,
                         pop_totals = NULL
                         )
                         {

  method <- get_method(method_selection)
  inv_link <- method$make_link_inv

  function(par) {
    theta <- as.matrix(par)
    X0 <- as.matrix(X)
    p <- ncol(X0)
    eta_pi <- X0 %*% theta
    ps <- inv_link(eta_pi)
    ps <- as.vector(ps)
    R_rand <- 1 - R

    R <- as.vector(R) # <------ required if Rcpp
    weights <- as.vector(weights) # <------ required if Rcpp
    R_rand <- as.vector(R_rand) # <------ required if Rcpp

    if (method_selection == "probit") {
      dinv_link <- method$make_link_inv_der
      psd <- dinv_link(eta_pi)
      psd <- as.vector(psd)
    }
    N_nons <- sum(1/ps)


    if (h == "1" || !is.null(pop_totals)) {
      mxDer <-  switch(method_selection,
                       "logit" = t(R * as.data.frame(X0) * (1-ps)/ps) %*% X0,
                       "cloglog" = t(R * as.data.frame(X0) * (1-ps)/ps^2 * exp(eta_pi)) %*% X0,
                       "probit" = t(R * as.data.frame(X0) * psd/ps^2) %*% X0)
    } else if (h == "2") {
      mxDer <- switch(method_selection,
                      "logit" = t(R_rand * as.data.frame(X0) * weights * ps/(exp(eta_pi)+1)) %*% X0,
                      "cloglog" = t(R_rand * as.data.frame(X0) * weights * (1-ps) * exp(eta_pi)) %*% X0,
                      "probit" = t(R_rand * as.data.frame(X0) * weights * psd) %*% X0)
    }
    matrix(mxDer/N_nons, nrow = p)

  }
}

q_lambda <- function(par,
                     lambda,
                     a = 3.7) {
  # SCAD penalty derivative
  penaltyd <- (abs(par)<lambda) * lambda + (abs(par)>=lambda) * ((a * lambda) > abs(par)) * ((a * lambda) - abs(par))/(a-1)
  penaltyd[1]<-0 # no penalty on the intercept
  penaltyd
}

# loss function for theta using the square distance of X between probability and nonprobability sample
# for lambda_theta selection

loss_theta <- function(par,
                       R,
                       X,
                       weights,
                       method_selection,
                       h) {

  method <- get_method(method_selection)

  inv_link <- method$make_link_inv

  theta <- as.matrix(par)
  nAB <- length(R)
  X0 <- as.matrix(X)
  eta_pi <- X0 %*% theta
  ps <- inv_link(eta_pi)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  R_rand <- 1 - R
  ps <- as.vector(ps)

  R <- as.vector(R) # <------ required if Rcpp
  weights <- as.vector(weights) # <------ required if Rcpp
  R_rand <- as.vector(R_rand) # <------ required if Rcpp

  loss <- switch(h,
                 "1" = sum(apply((X0*R/ps - X0*R_rand*weights), 2, sum)^2),
                 "2" = sum(apply((X0*R - X0*R_rand*weights*ps), 2, sum)^2))
  loss
}
