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
#' @importFrom ncvreg cv.ncvreg
#' @importFrom rootSolve multiroot
#' @importFrom stats qnorm
#' @importFrom survey svymean
#' @importFrom Matrix Matrix
#' @importFrom stats terms
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

  method <- get_method(method_selection)

  inv_link <- method$make_link_inv

  y_rand <- vector(mode = "numeric", length = n_rand)
  y <- c(y_rand, y_nons) # outcome variable for joint model
  n <- nrow(X)
  p <- ncol(X)
  N_rand <- sum(weights_rand)

  eps <- control_selection$epsilon
  maxit <- control_selection$maxit
  h <- control_selection$h_x
  lambda <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  nfolds <- control_selection$nfolds
  cv <- cv_nonprobsvy(X = X_stand, R = R, weights_X = weights_X,
                      method_selection = method_selection,
                      h = h, maxit = maxit, eps = eps,
                      lambda_min = lambda_min, nlambda = nlambda, nfolds = nfolds,
                      lambda = lambda)
  theta_est <- cv$theta_est[cv$theta_est != 0]
  min <- cv$min
  lambda <- cv$lambda
  theta_selected <- cv$theta_selected
  names(theta_est) <- c("(Intercept)", nons_names[theta_selected[-1]])


  # variables selection for beta using ncvreg package, perhaps cv.glmnet for theta like

  #theta <- glmnet::cv.glmnet(X, R, family = "binomial")
  #print(coef(theta))

  nlambda <- control_selection$nlambda
  penalty <- control_outcome$penalty
  beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1], y = y_nons,
                            penalty = penalty, family = family_outcome,
                            nlambda = nlambda)


  beta_est <- beta$fit$beta[,beta$min]
  beta_selected <- which(abs(beta_est) != 0) - 1
  beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

  # Estimating theta, beta parameters using selected variables

  idx <- unique(c(beta_selected[-1], theta_selected[-1])) #excluding intercepts
  psel <- length(idx)
  Xsel <- as.matrix(X[, idx + 1])

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
    svydesign_mean <- survey::svymean(~y_rand, svydesign)

    var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
    var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand * sigmasqhat))/N_nons^2 # nonprobability component

    se_nonprob <- sqrt(var_nonprob)
    se_prob <- sqrt(var_prob)

    var <- var_prob + var_nonprob #variance of an estimator

    se <- sqrt(var) # standard error

  } else if (family_outcome == "binomial") {

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
    svydesign_mean <- survey::svymean(~pi_rand, svydesign) # probability componemt based on survey package

    var_prob <- as.vector(attr(svydesign_mean, "var"))
    var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand*sigmasqhat))/(N_nons^2)

    se_nonprob <- sqrt(var_nonprob)
    se_prob <- sqrt(var_prob)


    # variance of an estimator
    var <- var_prob + var_nonprob

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

    svydesign_mean <- survey::svymean(~y_rand, svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var"))

    infl1 <- (y - y_hat)^2 * R/(ps_nons_est^2)
    infl2 <- (y - y_hat)^2 * R/ps_nons_est

    var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand*sigmasqhat))/(N_nons^2)

    se_nonprob <- sqrt(V2)
    se_prob <- sqrt(V1_svy)

    # variance of an estimator
    var <- var_prob + var_nonprob

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
         theta_sel = theta_est,
         beta_sel = beta_est,
         theta = theta_sel,
         beta = beta_sel
         ),
    class = "Doubly-robust")


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

  p <- ncol(X)
  n0 <- length(R)
  theta <- par[1:(p+1)]
  beta <- par[(p+2):(2*p+2)]
  X0 <- cbind(1, X)
  eta_pi <- X0 %*% theta
  ps <- inv_link(eta_pi)
  R_rand <- 1 - R
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  if (method_selection == "probit") {
    dinv_link <- method$make_link_inv_der
    psd <- dinv_link(eta_pi)
    psd <- as.vector(psd)
  }

  if (family_outcome == "gaussian") {

    res <- (y - (X0 %*% beta))
    res <- as.vector(res)

    UTB <- switch(method_selection,
                  "logit" = c(apply(X0*R/ps-X0*R_rand*weights, 2, sum), # estimating function
                              apply(X0*R*(1/ps-1)*res, 2, sum))/n0,
                  "cloglog" = c(apply(X0*R/ps-X0*R_rand*weights, 2, sum), # estimating function
                                apply(X0*R*(1-ps)/ps^2 * as.vector(exp(eta_pi))*res, 2, sum))/n0,
                  "probit" = c(apply(X0*R/ps-X0*R_rand*weights, 2, sum), # estimating function
                               apply(X0*R*psd/ps^2*res, 2, sum))/n0)

  } else if (family_outcome == "binomial") {

    m <- as.vector(exp(X0 %*% beta)/(1 + exp(X0 %*% beta)))
    res <- as.vector(y - m)
    m_der <- m * (1 - m) #derivative of m

    UTB <- switch(method_selection,
                  "logit" = c(apply(X0*R/ps*m_der - X0*R_rand*weights*m_der, 2, sum),
                              apply(X0*R*(1/ps-1)*res, 2, sum))/n0,
                  "cloglog" = c(apply(X0*R/ps*m_der - X0*R_rand*weights*m_der, 2, sum),
                                apply(X0*R*(1-ps)/ps^2 * exp(eta_pi)*res, 2, sum))/n0,
                  "probit" = c(apply(X0*R/ps*m_der - X0*R_rand*weights*m_der, 2, sum),
                               apply(X0*R*psd/ps^2*res, 2, sum))/n0)


  } else if (family_outcome == "poisson") {

    m <- as.vector(exp(X0 %*% beta))
    res <- as.vector(y - m)
    #derivative of is equal to m

    UTB <- switch(method_selection,
                "logit" = c(apply(X0*R/ps*m - X0*R_rand*weights*m, 2, sum),
                            apply(X0*R*(1/ps-1)*res, 2, sum))/n0,
                "cloglog" = c(apply(X0*R/ps*m - X0*R_rand*weights*m, 2, sum),
                              apply(X0*R*(1-ps)/ps^2 * exp(eta_pi)*res, 2, sum))/n0,
                "probit" = c(apply(X0*R/ps*m - X0*R_rand*weights*m, 2, sum),
                             apply(X0*R*psd/ps^2*res, 2, sum))/n0)

  }

  UTB

}
