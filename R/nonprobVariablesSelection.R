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

  lambda_theta <- control_selection$lambda
  lambda_beta <- control_outcome$lambda

  eps <- control_selection$epsilon

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample
  X_rand <- model.matrix(selection, svydesign$variables) #matrix of probability sample
  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  r_nons <- rep(1, nrow(X_nons))
  r_rand <- rep(0, nrow(X_rand))
  R_nons <- c(r_nons, r_rand)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise


  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand) # joint matrix
  y <- c(y_nons, rep(NA, n_rand)) # outcome variable for joint model
  y[which(is.na(y))] <- 0
  loc_nons <- which(R_nons == 1)
  loc_rand <- which(R_nons == 0)
  n <- dim(X)[1]
  p <- dim(X)[2]
  sw <- c(weights, d_rand) # vector of weights
  N_rand <- sum(d_rand)


  # initial values for set of parameters
  init_theta <- rep(0, p+1)
  init_beta <- rep(0, p+1)

  # variables selection using score equation for theta

  par0 <- c(init_theta)
  LAMBDA <- matrix(0,p+1,p+1)
  it <- 0

  for(jj in 1:100) {
    it <- it + 1

    Utheta0 <- Utheta(par = par0, R = R_nons, X = X, y = y,
                      d = d_rand, weights = weights,
                      method_selection = method_selection)

    Utheta0_der <- UthetaDer(par = par0, R = R_nons, X = X, y = y,
                             d = d_rand, weights = weights,
                             method_selection = method_selection)

    diag(LAMBDA) <- abs(q_lambda(par0, lambda_theta))/(eps + abs(par0))
    par <- par0 + MASS::ginv(Utheta0 + LAMBDA) %*% (Utheta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'

    if (sum(abs(par - par0)) < eps) break;
    if (sum(abs(par - par0)) > 1000) break;

    par0 <- par

  }

  par[which(abs(par) < 2 * eps)] <- 0
  theta_est <- par
  theta_selected <- which(theta_est != 0)


    # variables selection for beta using nvreg package

    beta <- ncvreg::ncvreg(X = X_nons, y = y_nons, lambda = lambda_beta,
                           penalty = 'SCAD', family = family_outcome)
    beta_est <- beta$beta
    beta_selected <- as.numeric(which(beta_est!=0))


    # Estimating theta, beta parameters using selected variables


    idx <- unique(c(beta_selected[-1] - 1, theta_selected[-1] - 1))
    psel <- length(idx)
    Xsel <- X[, idx]

    par0 <- rep(0, 2*(psel+1))

    # root for joint score equation
    par_sel <- rootSolve::multiroot(UThetaBeta,
                                    start = par0,
                                    R = R_nons,
                                    X = Xsel,
                                    y = y,
                                    d = d_rand,
                                    weights = weights,
                                    method_selection = method_selection,
                                    family_outcome = family_outcome)$root


    theta_sel <- par_sel[1:(psel+1)]
    beta_sel <- par_sel[(psel+2):(2*psel+2)]

    ps_nons_est  <- inv_link(as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(theta_sel)))
    d_nons <- 1/ps_nons_est[loc_nons]
    N_nons <- sum(d_nons)
    sw_rand <- sw[loc_rand]

    if (family_outcome == "gaussian") {

      y_hat <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = y_hat[loc_nons],
                         y_rand = y_hat[loc_rand],
                         d_nons = d_nons,
                         d_Rand = d_rand,
                         N_nons = N_nons,
                         N_rand = N_rand) #DR estimator

      # mu_hatdr <- sum((y - y_hat)*Rnons/ps_nons_est)/Nnons + (sum(y_hat*(1-Rnons) * sw.rand))/Nrand  # using mu_hatDR function in near future
      # using Nnons, Nrand instead of N

      y_rand <- y_hat[loc_rand]


      # sigmasqhat <- mean((y[loc_nons] - y_hat[loc_nons])^2) # squared errors mean
      # V1 <- sum((sw_rand^2 - sw_rand) * (y_rand)^2)/N_nons^2
      # V2 <- (sum(R_nons*(1-2*ps_nons_est)/ps_nons_est^2) + sum(sw_rand)) * sigmasqhat/(N_nons^2) #sum(sw_rand) = N_nons


      sigmasqhat <- mean((y[loc_rand] - y_hat[loc_rand])^2) # squared errors mean

      infl1 <- (y - y_hat)^2 * R_nons/(ps_nons_est^2)
      infl2 <- (y - y_hat)^2 * R_nons/ps_nons_est

      V1 <- sum((sw_rand^2 - sw_rand) * (y_rand)^2)/N_nons^2
      V2 <- (sum((infl1) - 2*infl2) + sum(sw_rand * sigmasqhat))/(N_nons^2)

      var <- V1 + V2 #variance of an estimator

      se <- sqrt(var) # standard error

    } else if (family_outcome == "binomial") {

      lm <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))
      pi <- exp(lm)/(1 + exp(lm))

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = pi[loc_nons],
                         y_rand = pi[loc_rand],
                         d_nons = d_nons,
                         d_rand = d_rand,
                         N_nons = N_nons,
                         N_rand = N_rand)

      # mu_hatdr <- sum((y - pi)*Rnons/ps_nons_est)/N + (sum(pi*(1-Rnons) * sw))/N  # using mu_hatDR function in near future
      # using Nnons, Nrand instead of N

      pi_rand <- pi[loc_rand]

      sigmasqhat <- pi[loc_rand] * (1 - pi[loc_rand])

      infl1 <- (y - pi)^2 * R_nons/(ps_nons_est^2)
      infl2 <- (y - pi)^2 * R_nons/ps_nons_est

      V1 <- sum((sw_rand^2 - sw_rand) * (pi_rand)^2)/N_nons^2
      V2 <- (sum((infl1) - 2*infl2) + sum(sw_rand*sigmasqhat))/(N_nons^2)


      # variance of an estimator
      var <- V1 + V2

      # standard error
      se <- sqrt(var)

    } else if (family_outcome == "poisson") { # to fix [variance]

      lm <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))
      y_hat <- exp(lm)

      mu_hat <- mu_hatDR(y = y_nons,
                         y_nons = y_hat[loc_nons],
                         y_rand = y_hat[loc_rand],
                         d_nons = d_nons,
                         d_rand = d_rand,
                         N_nons = N_nons,
                         N_rand = N_rand) #DR estimator

      # mu_hatdr <- sum((y - y_hat)*Rnons/ps_nons_est)/N + (sum(y_hat*(1-Rnons) * sw))/N

      sigmasqhat <- mean(y_hat[loc_rand])
      y_rand <- y_hat[loc_rand]

      infl1 <- (y - y_hat)^2 * R_nons/(ps_nons_est^2)
      infl2 <- (y - y_hat)^2 * R_nons/ps_nons_est

      V1 <- sum((sw_rand^2 - sw_rand) * (y_rand)^2)/N_nons^2
      V2 <- (sum((infl1) - 2*infl2) + sum(sw_rand*sigmasqhat))/(N_nons^2)

      # variance of an estimator
      var <- V1 + V2

      # standard error
      se <- sqrt(var)

    }

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    ci <- c(mu_hat - z * se, mu_hat + z * se)


  structure(
    list(population_mean = mu_hat,
         variance = var,
         CI = ci,
         standard_error = se,
         theta = theta_sel,
         beta = beta_sel),
    class = "Doubly-robust")


}



# score equation for theta, used in variable selection

Utheta <- function(par,
                   R,
                   X,
                   y,
                   d,
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
    sw <- c(weights, d)

    eq <- c(apply(X0 * R/ps - X0 * R_rand * sw, 2, sum))/n

    eq

}


# derivative of score equation for theta, used in variable selection

UthetaDer <-  function(par,
                       R,
                       X,
                       y,
                       d,
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
  sw <- rbind(weights, d)
  sw <- as.vector(sw)

  mxDer <- matrix(0,(p+1),(p+1))

  for (ii in 1:(p+1)) {
    for (jj in ii:(p+1)) {
      mxDer[ii,jj] <- sum(R * (1-ps)/ps * X0[,ii] * X0[,jj])
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

UThetaBeta <- function(par,
                       R,
                       X,
                       y,
                       d,
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
  sw <- c(weights, d)

  if (family_outcome == "gaussian") {

    res <- (y - (X0 %*% beta))
    res <- as.vector(res)

    UTB <- c(apply(X0*R/ps-X0*R_rand*sw, 2, sum), # estimating function
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  } else if (family_outcome == "binomial") {

    m <- exp(X0 %*% beta)/(1 + exp(X0 %*% beta))
    res <- (y - m)
    m_der <- m * (1 - m) #derivative of m

    UTB <- c(apply(X0*R/ps*m_der - X0*R_rand*sw*m_der, 2, sum),
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  } else if (family_outcome == "poisson") {

    m <- exp(X0 %*% beta)
    res <- (y - m)
    #derivative of is equal to m

    UTB <- c(apply(X0*R/ps*m - X0*R_rand*sw*m, 2, sum),
             apply(X0*R*(1/ps-1)*res, 2, sum))/n0

  }

  UTB

}


# loss function for theta using the square distance of X between probability and nonprobability sample

loss_theta <- function(par,
                       R,
                       y,
                       d,
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
  nAB <- length(Rnons)
  X0 <- cbind(1,X0)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)

  R_rand <- 1 - R
  ps <- as.vector(ps)
  sw <- rbind(weights, d)
  sw<-as.vector(sw)
  N_est_rand <- sum(d)
  N_est_nons <- sum(1/ps)

  loss <- sum(apply((x0*R/ps/N_est_nons - X0*R_rand*sw/N_est_rand), 2, sum)^2)
  loss

}



