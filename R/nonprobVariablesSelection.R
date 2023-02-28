# Implementation is based on IntegrativeFPM package, see at https://github.com/shuyang1987/IntegrativeFPM

#' Title nonprobSel
#'
#'
#' nonprobSel: Function for selecting important variables for sampling score model and outcome model
#'
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals - an optional `named vector` with population totals.
#' @param pop_means - an optional `named vector` with population means.
#' @param pop_size - an optional `double` with population size.
#' @param method_selection - a `character` with method for propensity scores estimation
#' @param method_outcome - a `character` with method for response variable estimation
#' @param family_selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - a
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a
#' @param control_selection a
#' @param control_outcome a
#' @param control_inference a
#' @param start a
#' @param verbose a
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @importFrom MASS ginv
#' @importFrom ncvreg ncvreg
#' @importFrom rootSolve multiroot
#' @importFrom stats qnorm
#' @importFrom survey svymean
#' @export
#'


nonprobSel <- function(selection,
                       outcome,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_selection,
                       method_outcome,
                       family_selection = "binomial",
                       family_outcome = "gaussian",
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection = controlSel(),
                       control_outcome = controlOut(),
                       control_inference = controlInf(),
                       start,
                       verbose,
                       contrasts,
                       model,
                       x,
                       y,
                       ...) {

  eps <- control_selection$epsilon

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  nons_names <- colnames(X_nons[,-1])
  if (all(nons_names %in% colnames(svydesign$variables))) {

    X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept

  } else {

    stop("variable names in data and svydesign do not match")

  }

  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_rand, R_nons)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_rand[, -1], X_nons[, -1]) # joint matrix
  weights_X <- c(weights_rand, weights)

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv

  y_rand <- rep(NA, n_rand)
  y <- c(y_rand, y_nons) # outcome variable for joint model
  y[which(is.na(y))] <- 0
  n <- dim(X)[1]
  p <- dim(X)[2]
  N_rand <- sum(weights_rand)

  lambda_theta <- cv_nonprobsvy(X_nons = X_nons,
                                X_rand = X_rand,
                                R = R,
                                weights_X = weights_X,
                                y = y,
                                method_selection = method_selection)


  # initial values for set of parameters
  init_theta <- rep(0, p+1)
  init_beta <- rep(0, p+1)

  # variables selection using score equation for theta

  par0 <- c(init_theta)
  LAMBDA <- matrix(0,p+1,p+1)
  it <- 0

  for(jj in 1:100) { # to fix
    it <- it + 1

    u_theta0 <- u_theta(par = par0, R = R, X = X, y = y,
                        weights = weights_X,
                        method_selection = method_selection)

    u_theta0_der <- u_theta_der(par = par0, R = R, X = X, y = y,
                               weights = weights_X,
                               method_selection = method_selection)

    diag(LAMBDA) <- abs(q_lambda(par0, lambda_theta))/(1e-6 + abs(par0))
    par <- par0 + MASS::ginv(u_theta0_der + LAMBDA) %*% (u_theta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'
                                                                                     # equation (13) in article
    if (sum(abs(par - par0)) < 1e-6) break;
    if (sum(abs(par - par0)) > 1000) break;

    par0 <- par

  }

  par[which(abs(par) < 1e-3)] <- 0
  theta_est <- par
  theta_selected <- which(theta_est != 0)


    # variables selection for beta using nvreg package

    #beta <- ncpen::ncpen(y.vec = y_nons, x.mat = X_nons, lambda = lambda_beta,
                         #penalty = "scad", family = family_outcome)

    beta <- ncvreg::cv.ncvreg(X = X[loc_nons, ], y = y_nons,             # perhaps using ncpen package -> these same results
                              penalty = 'SCAD', family = family_outcome)    # in ncvfit no possibility to define family

    beta_est <- beta$fit$beta[,beta$min]
    beta_selected <- as.numeric(which(abs(beta_est) > 2 * eps))

    # Estimating theta, beta parameters using selected variables

    idx <- unique(c(beta_selected[-1] - 1, theta_selected[-1] - 1)) #excluding intercepts
    psel <- length(idx)
    Xsel <- as.matrix(X[, idx])
    print(idx)

    par0 <- rep(0, 2*(psel+1))

    # root for joint score equation
    par_sel <- rootSolve::multiroot(u_theta_beta,
                                    start = par0,
                                    R = R,
                                    X = Xsel,
                                    y = y,
                                    weights = weights_X,
                                    method_selection = method_selection,
                                    family_outcome = family_outcome)$root


    theta_sel <- par_sel[1:(psel+1)]
    beta_sel <- par_sel[(psel+2):(2*psel+2)]
    # problem with theta_sel, beta_sel

    ps_nons_est <- inv_link(as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(theta_sel)))
    weights_nons <- 1/ps_nons_est[loc_nons]
    N_nons <- sum(weights_nons)

    if (family_outcome == "gaussian") {

      y_hat <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = y_hat[loc_nons],
                         y_rand = y_hat[loc_rand],
                         weights_nons = weights_nons,
                         weights_rand = weights_rand,
                         N_nons = N_nons,
                         N_rand = N_rand) #DR estimator

      y_rand <- y_hat[loc_rand]
      sigmasqhat <- mean((y[loc_rand] - y_hat[loc_rand])^2)

      infl1 <- (y - y_hat)^2 * R/(ps_nons_est^2)
      infl2 <- (y - y_hat)^2 * R/ps_nons_est

      # Variance estimatiors ####
      svydesign <- stats::update(svydesign,
                                 y_rand = y_rand)

      svydesign_mean <- survey::svymean(~y_rand, svydesign) # perhaps using survey package to compute prob variance
      V1_svy <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component

      V2 <- (sum((infl1) - 2*infl2) + sum(weights_rand * sigmasqhat))/(N_nons^2) # nonprobability component

      se_nonprob <- sqrt(V2)
      se_prob <- sqrt(V1_svy)

      var <- V1_svy + V2 #variance of an estimator

      se <- sqrt(var) # standard error

    } else if (family_outcome == "binomial") { # to do with Hajek approximation and/or survey package

      lm <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))
      pi <- exp(lm)/(1 + exp(lm))

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = pi[loc_nons],
                         y_rand = pi[loc_rand],
                         weights_nons = weights_nons,
                         weights_rand = weights_rand,
                         N_nons = N_nons,
                         N_rand = N_rand)


      pi_rand <- pi[loc_rand]

      sigmasqhat <- pi[loc_rand] * (1 - pi[loc_rand])

      infl1 <- (y - pi)^2 * R/(ps_nons_est^2)
      infl2 <- (y - pi)^2 * R/ps_nons_est

      # Variance estimators

      svydesign_mean <- survey::svymean(~pi_rand, svydesign) # perhaps using survey package to compute prob variance
                                                             # probability componemt based on survey package
      V1_svy <- var_prob <- attr(svydesign_mean, "var")

      V2 <- (sum((infl1) - 2*infl2) + sum(weights_rand*sigmasqhat))/(N_nons^2)

      se_nonprob <- sqrt(V2)
      se_prob <- sqrt(V1_svy)


      # variance of an estimator
      var <- as.vector(V1_svy + V2)

      # standard error
      se <- sqrt(var)

    } else if (family_outcome == "poisson") { # to fix [variance] # to do with Hajek approximation and/or survey package

      lm <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))
      y_hat <- exp(lm)

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = y_hat[loc_nons],
                         y_rand = y_hat[loc_rand],
                         weights_nons = weights_nons,
                         weights_rand = weights_rand,
                         N_nons = N_nons,
                         N_rand = N_rand) #DR estimator


      sigmasqhat <- mean(y_hat[loc_rand])
      y_rand <- y_hat[loc_rand]

      # Variance estimators



      svydesign <- stats::update(svydesign,
                                 y_rand = y_rand)

      svydesign_mean <- survey::svymean(~y_rand, svydesign) # perhaps using survey package to compute prob variance
      V1_svy <- attr(svydesign_mean, "var")

      infl1 <- (y - y_hat)^2 * R/(ps_nons_est^2)
      infl2 <- (y - y_hat)^2 * R/ps_nons_est

      V2 <- (sum((infl1) - 2*infl2) + sum(weiights_rand*sigmasqhat))/(N_nons^2)

      se_nonprob <- sqrt(V2)
      se_prob <- sqrt(V1_svy)

      # variance of an estimator
      var <- as.vector(V1_svy + V2)

      # standard error
      se <- sqrt(var)

    }

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    ci <- c(mu_hat - z * se, mu_hat + z * se)


  structure(
    list(mean = mu_hat,
         #VAR = var,
         CI = ci,
         SE = se,
         SE_nonprob = se_nonprob,
         SE_prob = se_prob,
         theta_est = theta_est,
         beta_est = beta_est,
         theta = theta_sel,
         beta = beta_sel,
         theta_selected = theta_selected[-1] - 1,
         beta_selected = beta_selected[-1] - 1
         ),
    class = "Doubly-robust")


}



# score equation for theta, used in variable selection

u_theta <- function(par,
                    R,
                    X,
                    y,
                    weights,
                    method_selection) {


    method <- method_selection
    if (is.character(method)) {
      method <- get(method, mode = "function", envir = parent.frame())
    }
    if (is.function(method)) {
      method <- method()
    }

    inv_link <- method$make_link_inv

    theta <- par
    n <- length(R)
    X0 <- cbind(1, X)
    lpiB <- X0 %*% theta
    ps <- inv_link(lpiB)
    R_rand <- 1 - R
    ps <- as.vector(ps)

    eq <- c(apply(X0 * R/ps - X0 * R_rand * weights, 2, sum))/n

    eq

}


# derivative of score equation for theta, used in variable selection

u_theta_der <-  function(par,
                         R,
                         X,
                         y,
                         weights,
                         method_selection) {

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv

  theta <- par
  p <- dim(X)[2]
  nAB <- length(R)
  X0 <- cbind(1, X)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)
  ps <- as.vector(ps)

  mxDer <- matrix(0,(p+1),(p+1))

  for (ii in 1:(p+1)) {
    for (jj in ii:(p+1)) {
      mxDer[ii,jj] <- mxDer[jj,ii] <- sum(R * (1-ps)/ps * X0[,ii] * X0[,jj])
    }
  }

  mxDer/nAB
}

q_lambda <- function(par,
                     lambda,
                     a = 3.7) {
  # SCAD penalty derivative

  penaltyd <- (abs(par)<lambda) * lambda + (abs(par)>=lambda) * ((a * lambda) > abs(par)) * ((a * lambda) - abs(par))/(a-1)
  penaltyd[1]<-0 # no penalty on the intercept

  penaltyd
}


# joint score equation for theta and beta, used in estimation

u_theta_beta <- function(par,
                         R,
                         X,
                         y,
                         weights,
                         method_selection,
                         family_outcome) {

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv


  p <- dim(X)[2]
  n0 <- length(R)
  theta <- par[1:(p+1)]
  beta <- par[(p+2):(2*p+2)]
  X0 <- cbind(1, X)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)
  R_rand <- 1 - R
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  if (family_outcome == "gaussian") {

    res <- (y - (X0 %*% beta))
    res <- as.vector(res)

    UTB <- c(apply(X0*R/ps-X0*R_rand*weights, 2, sum), # estimating function
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  } else if (family_outcome == "binomial") {

    m <- as.vector(exp(X0 %*% beta)/(1 + exp(X0 %*% beta)))
    res <- as.vector(y - m)
    m_der <- m * (1 - m) #derivative of m

    UTB <- c(apply(X0*R/ps*m_der - X0*R_rand*weights*m_der, 2, sum),
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  } else if (family_outcome == "poisson") {

    m <- as.vector(exp(X0 %*% beta))
    res <- as.vector(y - m)
    #derivative of is equal to m

    UTB <- c(apply(X0*R/ps*m - X0*R_rand*weights*m, 2, sum),
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  }

  UTB

}


# loss function for theta using the square distance of X between probability and nonprobability sample
# for selecting lambda_theta

loss_theta <- function(par,
                       R,
                       X,
                       y,
                       weights,
                       method_selection) {

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv


  theta <- par
  nAB <- length(R)
  X0 <- cbind(1,X)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  R_rand <- 1 - R
  ps <- as.vector(ps)
  N_est_rand <- sum(weights[loc_rand])
  N_est_nons <- sum(1/ps)

  loss <- sum(apply((X0*R/ps - X0*R_rand*weights), 2, sum)^2)

  loss

}

cv_nonprobsvy <- function(X_rand, X_nons, R, weights_X, y, method_selection, K = 10) {

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  R_nons <- R[loc_nons]
  R_rand <- R[loc_rand]

  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]

  weights_nons <- weights_X[loc_nons]
  weights_rand <- weights_X[loc_rand]

  k <- 1

  loss_theta_av <- vector(mode = "numeric", length = 50)
  lambdas <- seq(from = 0.02, to = 1, by = 0.02)

  X_nons <- X_nons[sample(nrow(X_nons)), ]
  X_rand <- X_rand[sample(nrow(X_rand)), ]

  folds_nons <- cut(seq(1,nrow(X_nons)), breaks=K, labels=FALSE) #split nonprobabability sample into K parts
  folds_rand <- cut(seq(1,nrow(X_rand)), breaks=K, labels=FALSE) #split probabability sample into K parts

  # pair K subsets randomly
  sample_nons <- sample(1:K, K, replace = FALSE)
  sample_rand <- sample(1:K, K, replace = FALSE)

  for (lambda in lambdas) {

    loss_theta_vec <- vector(mode = "numeric", length = K)

    for(i in 1:K){

      # train data for X_nons
      idx_nons <- which(folds_nons==sample_nons[i], arr.ind=TRUE)
      X_nons_train <- X_nons[-idx_nons, ]
      R_nons_train <- R_nons[-idx_nons]
      weights_nons_train <- weights_nons[-idx_nons]
      y_nons_train <- y[-idx_nons]

      # test data for X_nons
      X_nons_test <- X_nons[idx_nons, ]
      R_nons_test <- R_nons[idx_nons]
      weights_nons_test <- weights_nons[idx_nons]
      y_nons_test <- y[idx_nons]

      # train data for X_rand
      idx_rand <- which(folds_rand==sample_rand[i], arr.ind=TRUE)
      X_rand_train <- X_rand[-idx_rand, ]
      R_rand_train <- R_rand[-idx_rand]
      weights_rand_train <- weights_rand[-idx_rand]
      y_rand_train <- y_rand[-idx_rand]

      # test data for X_rand
      X_rand_test <- X_rand[idx_rand, ]
      R_rand_test <- R_rand[idx_rand]
      weights_rand_test <- weights_rand[idx_rand]
      y_rand_test <- y_rand[idx_rand]


      X_train <- rbind(X_rand_train[, -1], X_nons_train[, -1])
      X_test <- rbind(X_rand_test[, -1], X_nons_test[, -1])

      R_train <- c(R_rand_train, R_nons_train)
      R_test <- c(R_rand_test, R_nons_test)

      y_train <- c(y_rand_train, y_nons_train)
      y_test <- c(y_rand_test, y_nons_test)

      weights_X_train <- c(weights_rand_train, weights_nons_train)
      weights_X_test <- c(weights_rand_test, weights_nons_test)

      p <- dim(X_test)[2]

      # initial values for set of parameters
      init_theta <- rep(0, p+1)

      # variables selection using score equation for theta

      par0 <- c(init_theta)
      LAMBDA <- matrix(0,p+1,p+1)
      it <- 0

      for(jj in 1:100) { # to fix
        it <- it + 1

        u_theta0 <- u_theta(par = par0, R = R_train, X = X_train, y = y_train,
                            weights = weights_X_train,
                            method_selection = method_selection)

        u_theta0_der <- u_theta_der(par = par0, R = R_train, X = X_train, y = y_train,
                                    weights = weights_X_train,
                                    method_selection = method_selection)

        diag(LAMBDA) <- abs(q_lambda(par0, lambda))/(1e-6 + abs(par0))
        par <- par0 + MASS::ginv(u_theta0_der + LAMBDA) %*% (u_theta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'
        # equation (13) in the article
        if (sum(abs(par - par0)) < 1e-6) break;
        if (sum(abs(par - par0)) > 1000) break;

        par0 <- par

      }

      par[which(abs(par) < 1e-3)] <- 0
      theta_est <- par
      theta_selected <- which(theta_est != 0)

      loss_theta_vec[i] <- loss_theta(par = theta_est, R = R_test, X = X_test, y = y_test,
                                      weights = weights_X_test, method_selection = method_selection)


    }

    loss_theta_av[k] <- mean(loss_theta_vec)

    k <- k + 1

  }

  i <- which.min(loss_theta_av)

  lambdas[i]

}




