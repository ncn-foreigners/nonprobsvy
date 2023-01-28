#' Title Variable Selection
#'
#' VariableSelection: Function for selecting important variables for sampling score model and outcome model
#'
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop.totals - an optional `named vector` with population totals.
#' @param pop.means - an optional `named vector` with population means.
#' @param pop.size - an optional `double` with population size.
#' @param method.selection - a `character` with method for propensity scores estimation
#' @param method.outcome - a `character` with method for response variable estimation
#' @param family.selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action a
#' @param control.selection a
#' @param control.outcome a
#' @param control.inference a
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
#' @export
#'


nonprobSel <- function(selection,
                       outcome,
                       data,
                       svydesign,
                       pop.totals,
                       pop.means,
                       pop.size,
                       method.selection,
                       method.outcome,
                       family.selection = "binomial",
                       family.outcome = "gaussian",
                       subset,
                       weights,
                       na.action,
                       control.selection = controlSel(),
                       control.outcome = controlOut(),
                       control.inference = controlInf(),
                       start,
                       verbose,
                       contrasts,
                       model,
                       x,
                       y,
                       ...){

  lambda_theta <- control.selection$lambda
  lambda_beta <- control.outcome$lambda

  eps <- control.selection$epsilon

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample
  X_rand <- model.matrix(selection, svydesign$variables) #matrix of probability sample
  y_nons = XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  rnons <- rep(1, nrow(X_nons))
  rrand <- rep(0, nrow(X_rand))
  Rnons <- c(rnons, rrand)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise


  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$linkInv

  n.nons <- nrow(X_nons)
  n.rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand) # joint matrix
  y <- c(y_nons, rep(NA, n.rand)) # outcome variable for joint model
  y[which(is.na(y))] <- 0
  loc.nons <- which(Rnons == 1)
  loc.rand <- which(Rnons == 0)
  n <- dim(X)[1]
  p <- dim(X)[2]
  sw <- c(weights, d_rand) # vector of weights
  N <- sum(d_rand)


  ## initial values for set of parameters
  init.theta <- rep(0, p+1)
  init.beta <- rep(0, p+1)

  ## variables selection using score equation for theta

  par0 <- c(init.theta)
  LAMBDA <- matrix(0,p+1,p+1)
  it <- 0

  for(jj in 1:100){
    it <- it + 1

    Utheta0 <- Utheta(par0, Rnons, X, y, d_rand, weights, method.selection)
    Utheta0der <- UthetaDer(par0, Rnons, X, y, d_rand, weights, method.selection)

    diag(LAMBDA) <- abs(q_lambda(par0, lambda_theta))/(eps + abs(par0))
    par <- par0 + MASS::ginv(Utheta0 + LAMBDA) %*% (Utheta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'

    if(sum(abs(par - par0)) < eps) break;
    if(sum(abs(par - par0)) > 1000) break;
    par0 <- par
  }

  par[which(abs(par) < 2 * eps)] <- 0
  theta.est <- par
  theta.selected <- which(theta.est != 0)

  #for family == gaussian in outcome family

  # variabless selection for beta using nvreg package

  beta <- ncvreg::ncvfit(X_nons, y_nons, lambda = lambda_beta, penalty = 'SCAD')
  beta.est <- beta$beta
  beta.selected <- as.numeric(which(beta.est!=0))


  # Estimating parameters theta, beta using selected variables


  idx <- unique(c(beta.selected[-1] - 1, theta.selected[-1] - 1))
  psel <- length(idx)
  Xsel <- X[, idx]

  par0 <- rep(0, 2*(psel+1))

  ## root for joint score equation
  par_sel <- rootSolve::multiroot(UThetaBeta,
                                 par0,
                                 Rnons=Rnons,
                                 X=Xsel,
                                 y=y,
                                 d=d_rand,
                                 weights = weights,
                                 method.selection = method.selection)$root


  theta_sel <- par_sel[1:(psel+1)]
  beta_sel <- par_sel[(psel+2):(2*psel+2)]

  ps_nons_est  <- inv_link(as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(theta_sel)))
  y_hat <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))

  mu_hatdr <- sum((y - y_hat)*Rnons/ps_nons_est)/N + (sum(y_hat*(1-Rnons) * sw))/N  # using mu_hatDR function in near future
                                                                                    # using Nnons, Nrand instead of N

  sw.rand <- sw[loc.rand]
  y_rand <- y_hat[loc.rand]

  sigmasqhat <- mean((y[loc.nons] - y_hat[loc.nons])^2) # errors mean

  ve.pdr <- sum((sw.rand^2 - sw.rand) * (y_rand)^2)/N^2 + (sum(Rnons*(1-2*ps_nons_est)/ps_nons_est^2) + N) * sigmasqhat/(N^2) #variance of error term
  se.pdr <- sqrt(ve.pdr) # standard error

  return(list(theta = theta_sel,
         beta = beta_sel,
         mu_hat = mu_hatdr,
         ve = ve.pdr,
         se = se.pdr))


}



### score equation for theta, using in variable selection

Utheta <- function(par,
                   Rnons,
                   X,
                   y,
                   d,
                   weights,
                   method.selection){


    method <- method.selection
    if (is.character(method)) {
      method <- get(method, mode = "function", envir = parent.frame())
    }
    if (is.function(method)) {
      method <- method()
    }

    inv_link <- method$linkInv

    theta <- par
    n <- length(Rnons)
    X0 <- cbind(1, X)
    lpiB <- X0 %*% theta
    ps <- inv_link(lpiB)
    Rrand <- 1 - Rnons
    ps <- as.vector(ps)
    sw <- c(weights, d)

    eq <- c(apply(X0 * Rnons/ps - X0 * Rrand * sw, 2, sum))/n
    return(eq)

}


## derivative of score equation for theta, using in variable selection

UthetaDer <-  function(par,
                       Rnons,
                       X,
                       y,
                       d,
                       weights,
                       method.selection){

  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$linkInv

  theta <- par
  p <- dim(X)[2]
  nAB <- length(Rnons)
  X0 <- cbind(1, X)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)
  Rrand <- 1 - Rnons
  ps <- as.vector(ps)
  sw <- rbind(weights, d)
  sw <- as.vector(sw)

  mxDer <- matrix(0,(p+1),(p+1))

  for(ii in 1:(p+1) ){
    for(jj in ii:(p+1) ){
      mxDer[ii,jj] <- sum(Rnons * (1-ps)/ps * X0[,ii] * X0[,jj])
    }
  }

  return(mxDer/nAB)
}

q_lambda <- function(par,
                     lambda,
                     a = 3.7){
  ## SCAD penalty derivative

  penaltyd <- (abs(par)<lambda) * lambda + (abs(par)>=lambda) * ((a * lambda) > abs(par)) * ((a * lambda) - abs(par))/(a-1)
  penaltyd[1]<-0 # no penalty on the intercept
  return(penaltyd)
}


## joint score equation for theta and beta, using in estimation

UThetaBeta <- function(par,
                       Rnons,
                       X,
                       y,
                       d,
                       weights,
                       method.selection){

  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$linkInv


  p <- dim(X)[2]
  n0 <- length(Rnons)
  theta <- par[1:(p+1)]
  beta <- par[(p+2):(2*p+2)]
  X0 <- cbind(1, X)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)
  Rrand <- 1 - Rnons
  y[which(is.na(y))] <- 0
  res <- (y - (X0 %*% beta))
  ps <- as.vector(ps)
  sw <- c(weights, d)
  res <- as.vector(res)

  UTB <- c(apply(X0*Rnons/ps-X0*Rrand*sw,2,sum), # estimating function
           apply(X0*Rnons*(1/ps-1)*res,2,sum))/n0

  return(UTB)

}


## loss function for theta using the square distance of X between probability and nonprobability sample

loss_theta <- function(par,
                       Rnons,
                       y,
                       d,
                       weights,
                       method.selection){

  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$linkInv


  theta <- par
  nAB <- length(Rnons)
  X0 <- cbind(1,X0)
  lpiB <- X0 %*% theta
  ps <- inv_link(lpiB)

  Rrand <- 1 - Rnons
  ps <- as.vector(ps)
  sw <- rbind(weights, d)
  sw<-as.vector(sw)
  Nest_rand <- sum(d)
  Nest_nons <- sum(1/ps)

  loss <- sum(apply((x0*Rnons/ps/Nest_nons - X0*Rrand*sw/Nest_rand), 2, sum)^2)
  return(loss)

}



