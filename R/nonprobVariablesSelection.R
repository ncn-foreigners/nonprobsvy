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
#' @importFrom nleqslv nleqslv
#' @importFrom stats qnorm
#' @importFrom stats delete.response
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
    xx <- paste("~", paste(nons_names, collapse = "+"))
    outcome_rand <- as.formula(paste(outcome[2], xx))
    X_rand <- model.matrix(delete.response(terms(outcome_rand)), svydesign$variables[, nons_names]) # bug if formula is y~. #matrix of probability sample with intercept X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #
  } else {
    stop("variable names in data and svydesign do not match")
  } # TODO with dot_check for factor variables

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
  prior_weights <- c(weights_rand, weights)

  method <- get_method(method_selection)
  inv_link <- method$make_link_inv

  y_rand <- vector(mode = "numeric", length = n_rand)
  y <- c(y_rand, y_nons) # outcome variable for joint model
  n <- nrow(X)
  p <- ncol(X)
  N_rand <- sum(weights_rand)

  # Cross-validation for variable selection
  cv <- cv_nonprobsvy_rcpp(X = X_stand, # TODO add weights
                           R = R,
                           weights_X = prior_weights,
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
  X_design <- cbind(1, Xsel)
  #start <- start_fit(Xsel,
  #                   R,
  #                   weights,
  #                   weights_rand,
  #                   method_selection)
  #par0 <- c(rep(0, psel + 2), start)
  par0 <- rep(0, 2*(psel + 1))

  # root for joint score equation
  #multiroot <- rootSolve::multiroot(u_theta_beta_dr,
  #                                  start = par0,
  #                                  R = R,
  #                                  X = Xsel,
  #                                  y = y,
  #                                  weights = prior_weights,
  #                                  method_selection = method_selection,
  #                                  family_nonprobsvy = family_nonprobsvy)
  #par_sel <- multiroot$root

  multiroot <- nleqslv::nleqslv(x = par0, # TODO add user-specified parameters to control functions
                                fn = u_theta_beta_dr,
                                method = "Newton", # TODO consider the method
                                global = "qline",
                                xscalm = "fixed",
                                jacobian = TRUE,
                                control = list(scalex = rep(1, length(par0))), # TODO algorithm did not converge in maxit iterations for cloglog
                                R = R,
                                X = Xsel,
                                y = y,
                                weights = prior_weights,
                                method_selection = method_selection,
                                family_nonprobsvy = family_nonprobsvy)
  par_sel <- multiroot$x
  if (multiroot$termcd %in% c(2:7, -10)) {
    switch(as.character(multiroot$termcd),
           "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv, but user must check if function values are acceptably small."),
           "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable new point when fitting selection model by nleqslv."),
           "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
           "5" = warning("ill-conditioned Jacobian when fitting selection model by nleqslv."),
           "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
           "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
           "-10" = warning("user specified Jacobian is incorrect when fitting selection model by nleqslv."))
  }

  theta_hat <- par_sel[1:(psel+1)]
  beta_hat <- par_sel[(psel+2):(2*psel+2)]
  names(theta_hat) <- names(beta_hat) <- c("(Intercept)", colnames(X_nons)[idx+1]) # nons_names[idx]
  df_residual <- nrow(Xsel) - length(theta_hat)

  ps <- inv_link(as.vector(X_design %*% as.matrix(theta_hat)))
  weights_nons <- 1/ps[loc_nons]
  N_nons <- sum(weights * weights_nons)

  # variance-covariance matrix for selection model
  V <- diag(ps * (1 - ps))
  vcov_selection <- solve(t(X_design) %*% V %*% X_design)
  theta_errors <- sqrt(diag(vcov_selection))


  eta <- as.vector(X_design %*% as.matrix(beta_hat))
  y_hat <- family_nonprobsvy$mu(eta)
  y_rand_pred <- y_hat[loc_rand]
  y_nons_pred <- y_hat[loc_nons]
  sigma <- family_nonprobsvy$variance(mu = y_rand_pred, y = y[loc_rand])
  residuals <- family_nonprobsvy$residuals(mu = y_rand_pred, y = y[loc_rand])

  # variance-covariance matrix for outcome model
  vcov_outcome <- sigma * solve(t(X_design) %*% X_design)
  beta_errors <- sqrt(diag(vcov_outcome))

  # TODO std for coefficients

  selection <- list(theta_hat = theta_hat,
                    grad = multiroot$fvec[1:(psel+1)], # TODO
                    hess = NULL, # TODO
                    ps_nons = ps,
                    variance_covariance = vcov_selection,
                    df_residual = df_residual,
                    log_likelihood = "NULL")

  outcome <- list(beta_hat = beta_hat,
                  grad = multiroot$f.root[(psel+2):(2*psel+2)],
                  hess = NULL, # TODO
                  variance_covariance = vcov_outcome,
                  df_residual = df_residual,
                  family = list(mu = y_hat,
                                var = sigma,
                                residuals = residuals),
                  log_likelihood = "NULL")

  mu_hat <- mu_hatDR(y = y_nons,
                     y_nons = y_nons_pred,
                     y_rand = y_rand_pred,
                     weights = weights,
                     weights_nons = weights_nons,
                     weights_rand = weights_rand,
                     N_nons = N_nons,
                     N_rand = N_rand) #DR estimator

  infl1 <- (prior_weights * (y - y_hat))^2 * R/(ps^2) # TODO add weights
  infl2 <- (prior_weights * (y - y_hat))^2 * R/ps # TODO add weights

  # Variance estimators ####
  svydesign <- stats::update(svydesign,
                             y_rand = y_rand_pred)
  svydesign_mean <- survey::svymean(~y_rand, svydesign)

  var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
  var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand * sigma))/N_nons^2 # nonprobability component
  se_prob <- sqrt(var_prob)
  se_nonprob <- sqrt(var_nonprob)


  SE_values <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))

  var <- var_prob + var_nonprob #variance of an estimator
  se <- sqrt(var) # standard error
  output <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))

  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  confidence_interval <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                              upper_bound = mu_hat + z * se))))

  parameters <- matrix(c(theta_hat, theta_errors),
                       ncol = 2,
                       dimnames = list(names(theta_hat),
                                       c("Estimate", "Std. Error")))

  beta <- matrix(c(beta_hat, beta_errors),
                 ncol = 2,
                 dimnames = list(names(theta_hat),
                                 c("Estimate", "Std. Error")))

  probabilities_summary <- summary(as.vector(ps[loc_nons]))
  prop_scores <- as.vector(ps)

  structure(
    list(X = X,
         prop_scores = prop_scores,
         weights = as.vector(weights_nons),
         control = list(control_selection = control_selection,
                        control_outcome = control_outcome,
                        control_inference = control_inference),
         output = output,
         confidence_interval = confidence_interval,
         SE_values = SE_values,
         parameters = parameters,
         beta = beta,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = N_nons,
         #df_residual = df_residual,
         outcome = outcome,
         selection = selection
         #log_likelihood = "NULL"
         ),
    class = c("nonprobsvy", "nonprobsvy_dr"))
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
  #weights <- rep.int(1, nrow(data)) # to remove

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

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
  y_rand <- vector(mode = "numeric", length = n_rand)
  y <- c(y_rand, y_nons) # outcome variable for joint model
  X <- rbind(X_rand, X_nons) # joint matrix
  X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
  prior_weights <- c(weights_rand, weights)

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
  X <- rbind(X_nons, X_rand) # joint model matrix

  # TODO Estimation with bias minimisation
  # multiroot <- rootSolve::multiroot(u_theta_beta_mi,
  #                                  start = rep(0, ncol(Xsel) + 1),
  #                                  R = R,
  #                                  X = Xsel,
  #                                  y = y,
  #                                  weights = weights_X,
  #                                  method_selection = method_selection,
  #                                  family_outcome = family_outcome)

  #beta_sel <- multiroot$root

  # Estimation for outcome model with ordinary least squares method (glm)
  model_out <- internal_outcome(outcome = outcome,
                                data = data[,beta_selected+1],
                                weights = weights,
                                family_outcome = family_outcome)

  beta_sel <- model_out$glm_summary$coefficients[,1]
  parameters <- model_out$glm_summary$coefficients
  names(beta_sel) <- c("(Intercept)", nons_names[beta_selected])
  #beta <- model_out$glm_summary$coefficients
  N_est_rand <- sum(weights_rand)
  y_rand_pred <- as.numeric(X_rand %*% beta_sel) # y_hat for probability sample # consider predict function
  y_nons_pred <- model_out$glm$fitted.values #as.numeric(X_nons %*% model_nons_coefs)
  df_residual <- nrow(X) - length(beta_sel)

  outcome <- list(beta_hat = beta_sel,
                  grad = NULL,
                  hess = NULL,
                  variance_covariance = vcov(model_out$glm),
                  df_residual = df_residual,
                  log_likelihood = "NULL")

  #par0 <- rep(0, length(beta_sel))

  #multiroot <- nleqslv::nleqslv(x = par0, # TODO add user-specified parameters to control functions
  #                              fn = u_theta_beta_mi,
  #                              method = "Newton", # TODO consider the method
  #                              global = "qline",
  #                              xscalm = "fixed",
  #                              jacobian = TRUE,
  #                              R = R,
  #                              X = Xsel,
  #                              y = y,
  #                              weights = prior_weights,
  #                              family_nonprobsvy = family_nonprobsvy)
  #print(multiroot$x)

  mu_hat <- mu_hatMI(y = y_rand_pred,
                     weights_rand = weights_rand,
                     N = N_est_rand)
  if (is.null(pop_size)) pop_size <- N_est_rand # estimated pop_size

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
    SE_values <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))


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
    SE_values <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))

  }

  se <- sqrt(var)
  output <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))

  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  confidence_interval <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                              upper_bound = mu_hat + z * se))))


  structure(
    list(X = X,
         control = list(control_outcome = control_outcome,
                        control_inference = control_inference),
         output = output,
         confidence_interval = confidence_interval,
         SE_values = SE_values,
         parameters = parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         outcome = outcome
    ),
    class = c("nonprobsvy", "nonprobsvy_mi"))
}


nonprobSelP <- function(selection, # TODO
                        target,
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_selection,
                        method_outcome,
                        family_selection,
                        family_outcome,
                        subset,
                        strata,
                        weights,
                        na_action,
                        control_selection,
                        control_outcome,
                        control_inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...) { #TODO

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

  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  nons_names <- attr(terms(outcome, data = data), "term.labels")
  if (all(nons_names %in% colnames(svydesign$variables))) {
    X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #X_rand <- model.matrix(delete.response(terms(outcome)) , svydesign$variables) bug if formula is y~. #matrix of probability sample with intercept
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

  # Cross-validation for variable selection
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

  idx <- theta_selected[-1] # excluding intercepts
  psel <- length(idx)
  Xsel <- as.matrix(X[, idx + 1])
  X_design <- cbind(1, Xsel)

  par1 <- rep(0, length(theta_est))
  par0 <- start_fit(Xsel,
                    R,
                    weights,
                    weights_rand,
                    method_selection)

  #multiroot <- rootSolve::multiroot(u_theta_beta_ipw,
  #                                  start = c(0, par0),
  #                                  R = R,
  #                                  X = Xsel,
  #                                  y = y,
  #                                  weights = weights_X,
  #                                  method_selection = method_selection)

  multiroot <- nleqslv::nleqslv(x = c(0, par0), # TODO add user-specified parameters to control functions
                                fn = u_theta_beta_ipw,
                                method = "Newton", # TODO consider the method
                                global = "qline",
                                xscalm = "fixed",
                                jacobian = TRUE,
                                R = R,
                                X = Xsel,
                                y = y,
                                weights = weights_X,
                                method_selection = method_selection)

  theta_hat <- multiroot$x

  names(theta_hat) <- c("(Intercept)", nons_names[idx])
  df_residual <- nrow(Xsel) - length(theta_hat)

  ps <- inv_link(as.vector(as.matrix(X_design) %*% as.matrix(theta_hat)))
  weights_nons <- 1/ps[loc_nons]
  N_nons <- sum(weights, weights_nons)

  # variance-covariance matrix for selection model
  V <- diag(ps * (1 - ps))
  vcov_selection <- solve(t(X_design) %*% V %*% X_design)
  theta_errors <- sqrt(diag(vcov_selection))

  selection <- list(theta_hat = theta_hat,
                    grad = multiroot$fvec, # TODO
                    hess = NULL, # TODO
                    ps_nons = ps[loc_nons],
                    variance_covariance = vcov_selection,
                    df_residual = df_residual,
                    log_likelihood = "NULL")

  mu_hat <- mu_hatIPW(y = y_nons,
                      weights = weights,
                      weights_nons = weights_nons,
                      N = N_nons)
  SE_values <- data.frame(t(data.frame("SE" = c(prob = .5, nonprob = .5))))
  se <- 1 # TODO now just small number
  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)
    # confidence interval based on the normal approximation
  confidence_interval <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                              upper_bound = mu_hat + z * se
  ))))
  parameters <- matrix(c(theta_hat, theta_errors),
                       ncol = 2,
                       dimnames = list(names(theta_hat),
                                       c("Estimate", "Std. Error")))

  output <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))

  structure(
    list(X = X_design,
         prop_scores = ps,
         weights = weights_nons,
         control = list(control_selection = control_selection,
                        control_inference = control_inference),
         output = output,
         SE = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         selection = selection
    ),
    class = c("nonprobsvy", "nonprobsvy_ipw"))


}

# joint score equation for theta and beta, used in estimation when variable selections
u_theta_beta_dr <- function(par,
                            R,
                            X,
                            y,
                            weights,
                            method_selection,
                            family_nonprobsvy) {

  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev

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

  n <- length(R)
  R_rand <- 1 - R

  utb <- c(apply(X0 * R/ps * mu_der * weights - X0 * R_rand * weights * mu_der, 2, sum),
           apply(X0 * R * weights * as.vector(-inv_link_rev(eta_pi)) * res, 2, sum))/n

  utb
}


u_theta_beta_ipw <- function(par,
                             R,
                             X,
                             y,
                             weights,
                             method_selection) { # TODO

  method <- get_method(method_selection)
  inv_link_rev <- method$make_link_inv_rev
  inv_link <- method$make_link_inv

  p <- ncol(X)
  theta <- par
  X0 <- cbind(1, X)
  eta_pi <- X0 %*% theta
  y[which(is.na(y))] <- 0

  R_rand <- 1- R
  n <- length(R)

  UTB <- apply(X0 * (R * as.vector(inv_link_rev(eta_pi)) * y - y), 2, sum)/n

  UTB

}

u_theta_beta_mi <- function(par,
                            R,
                            X,
                            y,
                            weights,
                            family_nonprobsvy) { # TODO

  if(is.character(family_nonprobsvy)) {
    family_nonprobsvy <- paste(family_nonprobsvy, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }


  p <- ncol(X)
  beta <- par
  eta <- X %*% beta
  mu <- family_nonprobsvy$mu(eta)
  mu_der <- family_nonprobsvy$mu_der(mu)

  n <- length(R)
  R_rand <- 1 - R

  UTB <- apply(X * (R_rand * weights * as.vector(mu) - y), 2, sum)/n # TODO

  UTB
}

# TODO Jacobian of the estimating equations for dr method
u_theta_beta_dr_jacob <- function(par,
                                  R,
                                  X,
                                  y,
                                  weights,
                                  method_selection,
                                  family_nonprobsvy){
  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev
  dinv_link_rev <- method$make_link_inv_rev_de

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
  mu_der2 <- family_nonprobsvy$mu_der2(mu)
  res <- family_nonprobsvy$residuals(mu = mu, y = y)
  n <- length(R)
  R_rand <- 1 - R

  jac <- c(apply(- X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum),
           apply(X0 * R/ps * mu_der2 * weights - X0 * R_rand * weights * mu_der2, 2, sum),
           apply(X0 * R * weights * as.vector(dinv_link_rev(eta_pi)) * res * X0, 2, sum),
           apply(X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum))/n
  jac

}
