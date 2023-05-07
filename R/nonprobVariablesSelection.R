#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprobSel} fits model for inference based on non-probability surveys using various methods
#' with variable selection techniques and doubly robust approach. Implementation is based on [IntegrativeFPM package](https://github.com/shuyang1987/IntegrativeFPM).
#'
#' \loadmathjax
#'
#' @param selection `formula`, the selection (propensity) equation.
#' @param outcome `formula`, the outcome equation.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals an optional `named vector` with population totals.
#' @param pop_means an optional `named vector` with population means.
#' @param pop_size an optional `double` with population size.
#' @param method_selection a `character` with method for propensity scores estimation
#' @param method_outcome a `character` with method for response variable estimation
#' @param family_selection a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probability and non-probability samples, contains parameters such as estimation method or variance method
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation.
#' @param verbose verbose, numeric.
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... Additional, optional arguments.
#'
#' @useDynLib nonprobsvy
#' @importFrom MASS ginv
#' @importFrom ncvreg cv.ncvreg
#' @importFrom rootSolve multiroot
#' @importFrom stats qnorm
#' @importFrom survey svymean
#' @importFrom Matrix Matrix
#' @importFrom stats terms
#' @import RcppArmadillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
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

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  #if(is.function(family_outcome)) family_outcome <- family_outcome()

  eps <- control_selection$epsilon
  maxit <- control_selection$maxit
  h <- control_selection$h_x
  lambda <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  nfolds <- control_selection$nfolds
  #weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  nons_names <- attr(terms(outcome, data = data), "term.labels")
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
  X <- rbind(X_rand, X_nons) # joint matrix
  X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
  weights_X <- c(weights_rand, weights)

  method <- get_method(method_selection)
  inv_link <- method$make_link_inv

  y_rand <- vector(mode = "numeric", length = n_rand)
  y <- c(y_rand, y_nons) # outcome variable for joint model
  n <- nrow(X)
  p <- ncol(X)
  N_rand <- sum(weights_rand)

  cv <- cv_nonprobsvy_rcpp(X = X_stand,
                           R = R,
                           weights_X = weights_X,
                           method_selection = method_selection,
                           h = h,
                           maxit = maxit,
                           eps = eps,
                           lambda_min = lambda_min,
                           nlambda = nlambda,
                           nfolds = nfolds,
                           lambda = lambda)
  theta_est <- cv$theta_est[cv$theta_est != 0]
  min <- cv$min
  lambda <- cv$lambda
  theta_selected <- cv$theta_selected
  names(theta_est) <- c("(Intercept)", nons_names[theta_selected[-1]])

  nlambda <- control_outcome$nlambda
  penalty <- control_outcome$penalty
  beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1],
                            y = y_nons,
                            penalty = penalty,
                            family = family_outcome,
                            nlambda = nlambda)


  beta_est <- beta$fit$beta[,beta$min]
  beta_selected <- which(abs(beta_est) != 0) - 1
  beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

  # Estimating theta, beta parameters using selected variables

  idx <- unique(c(beta_selected[-1], theta_selected[-1])) # excluding intercepts
  psel <- length(idx)
  Xsel <- as.matrix(X[, idx + 1])
  #start <- start_fit(Xsel,
  #                   R,
  #                   weights,
  #                   weights_rand,
  #                   method_selection)
  #par0 <- c(rep(0, psel + 2), start)
  par0 <- rep(0, 2*(psel + 1))

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
  names(theta_sel) <- names(beta_sel) <- c("(Intercept)", nons_names[idx])

  ps_nons_est <- inv_link(as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(theta_sel)))
  weights_nons <- 1/ps_nons_est[loc_nons]
  N_nons <- sum(weights_nons)

  eta <- as.vector(as.matrix(cbind(1, Xsel)) %*% as.matrix(beta_sel))
  y_hat <- family_nonprobsvy$mu(eta)
  y_rand <- y_hat[loc_rand]
  y_nons <- y_hat[loc_nons]
  sigma <- family_nonprobsvy$variance(mu = y_rand, y = y[loc_rand])

  mu_hat <- mu_hatDR(y = y_nons,
                     y_nons = y_nons,
                     y_rand = y_rand,
                     weights_nons = weights_nons,
                     weights_rand = weights_rand,
                     N_nons = N_nons,
                     N_rand = N_rand) #DR estimator

  infl1 <- (y - y_hat)^2 * R/(ps_nons_est^2)
  infl2 <- (y - y_hat)^2 * R/ps_nons_est

  # Variance estimatiors ####
  svydesign <- stats::update(svydesign,
                             y_rand = y_rand)
  svydesign_mean <- survey::svymean(~y_rand, svydesign)

  var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
  var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand * sigma))/N_nons^2 # nonprobability component

  se_nonprob <- sqrt(var_nonprob)
  se_prob <- sqrt(var_prob)

  var <- var_prob + var_nonprob #variance of an estimator
  se <- sqrt(var) # standard error

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
         theta_sel = theta_est,
         beta_sel = beta_est,
         theta = theta_sel,
         beta = beta_sel
         ),
    class = "Doubly-robust")
}

#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprobSelM} fits model for inference based on non-probability surveys using various methods
#' with variable selection techniques and mass imputation approach.
#'
#' \loadmathjax
#'
#' @param outcome `formula`, the outcome equation.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals an optional `named vector` with population totals.
#' @param pop_means an optional `named vector` with population means.
#' @param pop_size an optional `double` with population size.
#' @param method_outcome a `character` with method for response variable estimation
#' @param family_outcome a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probability and non-probability samples, contains parameters such as estimation method or variance method
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation.
#' @param verbose verbose, numeric.
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... Additional, optional arguments.
#'


nonprobSelM <- function(outcome,
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_outcome,
                        family_outcome,
                        subset,
                        strata,
                        weights,
                        na_action,
                        control_outcome,
                        control_inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...) {
  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  nons_names <- attr(terms(outcome, data = data), "term.labels")
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
  X <- rbind(X_rand, X_nons) # joint matrix
  X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
  weights_X <- c(weights_rand, weights)

  nlambda <- control_outcome$nlambda
  penalty <- control_outcome$penalty
  beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1],
                            y = y_nons,
                            penalty = penalty,
                            family = family_outcome,
                            nlambda = nlambda)


  beta_est <- beta$fit$beta[,beta$min]
  beta_selected <- which(abs(beta_est) != 0) - 1
  beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

  Xsel <- as.matrix(X[, beta_selected + 1])

  X_rand <- Xsel[loc_rand, ]
  X_nons <- Xsel[loc_nons, ]

  # Estimation for outcome model
  model_out <- internal_outcome(X_nons = X_nons,
                                X_rand = X_rand,
                                y = y_nons,
                                weights = weights,
                                family_outcome = family_outcome)

  y_rand_pred <- model_out$y_rand_pred
  y_nons_pred <- model_out$y_nons_pred
  beta_sel <- model_out$model_nons_coefs
  names(beta_sel) <- c("(Intercept)", nons_names[beta_selected])
  N_est_rand <- sum(weights_rand)

  mu_hat <- mu_hatMI(y = y_rand_pred,
                     weights = weights_rand,
                     N = N_est_rand)

  if (control_inference$var_method == "analytic"){

    # updating probability sample by adding y_hat variable
    svydesign <- stats::update(svydesign,
                               y_hat_MI = y_rand_pred)

    svydesign_mean <- survey::svymean(~y_hat_MI, svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var"))

    mx <- 1/N_est_rand * colSums(weights_rand * X_rand)
    c <- solve(1/n_nons * t(X_nons) %*% X_nons) %*% mx
    e <- y_nons - y_nons_pred

    # nonprobability component
    var_nonprob <- 1/n_nons^2 * t(as.matrix(e^2)) %*% (X_nons %*% c)^2
    var_nonprob <- as.vector(var_nonprob)

    se_nonprob <- sqrt(var_nonprob)
    se_prob <- sqrt(var_prob)
    # variance
    var <- var_nonprob + var_prob

  } else if (control_inference$var_method == "bootstrap") {
    # bootstrap variance
    var <- bootMI(X_rand,
                  X_nons,
                  weights,
                  y_nons,
                  family_outcome,
                  1000,
                  weights_rand,
                  mu_hat,
                  svydesign,
                  rep_type = control_inference$rep_type)
    inf <- "not computed for bootstrap variance"

  }

  se <- sqrt(var)

  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  ci <- c(mu_hat - z*se, mu_hat + z*se)

  structure(
    list(mean = mu_hat,
         #VAR = var,
         CI = ci,
         SE = se,
         SE_nonprob = ifelse(control_inference$var_method == "analytic", se_nonprob, inf),
         SE_prob = ifelse(control_inference$var_method == "analytic", se_prob, inf),
         beta = beta_sel
    ),
    class = "Mass Imputation")
}

# joint score equation for theta and beta, used in estimation

u_theta_beta <- function(par,
                         R,
                         X,
                         y,
                         weights,
                         method_selection,
                         family_outcome) {

  method <- get_method(method_selection)

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  inv_link <- method$make_link_inv

  p <- ncol(X)
  theta <- par[1:(p+1)]
  beta <- par[(p+2):(2*p+2)]
  X0 <- cbind(1, X)
  eta_pi <- X0 %*% theta
  ps <- inv_link(eta_pi)
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  eta <- X0 %*% beta
  mu <- family_nonprobsvy$mu(eta)
  mu_der <- family_nonprobsvy$mu_der(mu)
  res <- family_nonprobsvy$residuals(mu = mu, y = y)


  UTB <- method$UTB(X = X0,
                    R = R,
                    weights = weights,
                    ps,
                    eta_pi = eta_pi,
                    mu_der = mu_der,
                    res = res)

  UTB
}
