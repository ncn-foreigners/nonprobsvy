# These functions are only used internally in the nonprobSel function, so there is no need for documenting them
#' @importFrom stats glm
#' @importFrom stats residuals

cv_nonprobsvy <- function(X, R, weights_X, method_selection, h, maxit, eps, lambda_min, nlambda, nfolds) {

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- cbind(X[loc_nons,], weights_X[loc_nons], R[loc_nons])
  X_rand <- cbind(X[loc_rand,], weights_X[loc_rand], R[loc_rand])
  k <- 1

  lambdas <- setup_lambda(X = X,
                          y = R,
                          weights = weights_X,
                          method_selection = method_selection,
                          lambda_min = 0,
                          nlambda = 50)
  #theta <- list(vector(mode = "numeric", length = length(lambdas)))
  loss_theta_av <- vector(mode = "numeric", length = length(lambdas))

  X_nons <- X_nons[sample(nrow(X_nons)), ]
  X_rand <- X_rand[sample(nrow(X_rand)), ]

  folds_nons <- cut(seq(1,nrow(X_nons)), breaks=nfolds, labels=FALSE) #split nonprobability sample into K parts
  folds_rand <- cut(seq(1,nrow(X_rand)), breaks=nfolds, labels=FALSE) #split probability sample into K parts

  # pair K subsets randomly
  sample_nons <- sample(1:nfolds, nfolds, replace = FALSE)
  sample_rand <- sample(1:nfolds, nfolds, replace = FALSE)

  for (lambda in lambdas) {

    loss_theta_vec <- vector(mode = "numeric", length = nfolds)

    for(i in 1:nfolds){
      # train data for X_nons
      idx_nons <- which(folds_nons==sample_nons[i], arr.ind=TRUE)
      X_nons_train <- X_nons[-idx_nons, ]
      # test data for X_nons
      X_nons_test <- X_nons[idx_nons, ]

      # train data for X_rand
      idx_rand <- which(folds_rand==sample_rand[i], arr.ind=TRUE)
      X_rand_train <- X_rand[-idx_rand, ]
      # test data for X_rand
      X_rand_test <- X_rand[idx_rand, ]

      X_train <- rbind(X_rand_train, X_nons_train)
      X_test <- rbind(X_rand_test, X_nons_test)
      ncols <- ncol(X_test)
      idxx <- 1:(ncols-2)

      theta_est <- fit_nonprobsvy(X_train[,idxx], R = X_train[, ncols],  weights = X_train[, ncols-1],
                                  method_selection, h, lambda, maxit, eps)
      if (any(theta_est==0)) {
        i <- which(theta_est==0)
        idxx <- idxx[-i]
      }
      X_testloss <- X_test[,idxx]
      R_testloss <- X_test[, ncols]
      weights_testloss <- X_test[, ncols-1]
      loss_theta_vec[i] <- loss_theta(par = theta_est[idxx], X = X_testloss, R = R_testloss,
                                      weights = weights_testloss, h = h, method_selection = method_selection)
    }
    loss_theta_av[k] <- mean(loss_theta_vec)
    k <- k + 1
  }
  min <- which.min(loss_theta_av)
  #print(data.frame(lambda = lambdas, loss = loss_theta_av))
  lambda_min <- lambdas[min]
  theta_est <- fit_nonprobsvy(X = X, R = R,  weights = weights_X,
                             method_selection, h, lambda_min, maxit, eps, warn = TRUE)
  theta_selected <- which(theta_est != 0) - 1
  list(min = min,
       lambda_min = lambda_min,
       theta_est = theta_est,
       theta_selected = theta_selected)

}



setup_lambda <- function(X, y, weights, method_selection, alpha = 1, lambda_min, log_lambda = FALSE, nlambda, ...) { #consider panalty factor here

  #fit <- glm.fit(x = X, y = y, weights = weights, family = binomial(link = method_selection))
  fit <- stats::glm(y~1, weights = weights, family = binomial(link = method_selection))

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


fit_nonprobsvy <- function(X, R, weights, method_selection, h, lambda, maxit, eps, warn = FALSE,  ...) {

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

    u_theta0 <- u_theta(R = R, X = X,
                        weights = weights, h = h,
                        method_selection = method_selection)

    u_theta0_der <- u_theta_der(R = R, X = X,
                                weights = weights, h = h,
                                method_selection = method_selection)

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
                    N = NULL) {


  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

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

    eq <- switch(h,
                 "1" = c(apply(X0 * R/ps - X0 * R_rand * weights, 2, sum))/N_nons,
                 "2" = c(apply(X0 * R - X0 * R_rand * weights * ps, 2, sum))/N_nons)
    eq
  }
}


# derivative of score equation for theta, used in variable selection

u_theta_der <-  function(R,
                         X,
                         weights,
                         method_selection,
                         h,
                         N = NULL) {

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv

  function(par) {
    theta <- as.matrix(par)
    X0 <- as.matrix(X)
    p <- ncol(X0)
    n <- length(R)
    eta_pi <- X0 %*% theta
    ps <- inv_link(eta_pi)
    ps <- as.vector(ps)
    R_rand <- 1 - R

    if (method_selection == "probit") {
      dinv_link <- method$make_link_inv_der
      psd <- dinv_link(eta_pi)
      psd <- as.vector(psd)
    }

    N_nons <- sum(1/ps)

    if (h == "1") {
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
    mxDer/N_nons

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
# for selecting lambda_theta

loss_theta <- function(par,
                       R,
                       X,
                       weights,
                       h,
                       method_selection) {

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

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
  N_est_rand <- sum(weights[loc_rand])
  N_est_nons <- sum(1/ps)

  loss <- switch(h,
                 "1" = sum(apply((X0*R/ps - X0*R_rand*weights), 2, sum)^2),
                 "2" = sum(apply((X0*R - X0*R_rand*weights*ps), 2, sum)^2))
  loss
}
