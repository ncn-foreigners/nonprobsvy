#' An internal function for bias correction Model for Non-probability Surveys
#'
#' @description
#' Fits a bias correction model for non-probability surveys using maximum likelihood estimation.
#'
#' @param X Matrix of covariates for both selection and outcome models
#' @param y Vector of outcome variables
#' @param weights Numeric vector of survey weights for non-probability sample
#' @param weights_rand Numeric vector of survey weights for probability sample
#' @param R Binary vector indicating sample membership (1 for non-probability, 0 for probability sample)
#' @param n_nons Integer, size of non-probability sample
#' @param n_rand Integer, size of probability sample
#' @param method_selection Character, specifies the selection model type
#' @param family Family object specifying the error distribution and link function
#' @param start_selection Numeric vector of starting values for selection model parameters
#' @param start_outcome Numeric vector of starting values for outcome model parameters
#' @param maxit Integer, maximum number of iterations for the optimization
#' @param nleqslv_method Character, method to be used in nleqslv optimization
#' @param nleqslv_global Character, global strategy for nleqslv
#' @param nleqslv_xscalm Character, scaling method for nleqslv
#' @param boot Logical, whether the function is being called for bootstrap (default FALSE)
#'
#' @return A list containing two components:
#' \itemize{
#'   \item selection: List of selection model results including:
#'     \itemize{
#'       \item theta_hat: Estimated selection model coefficients
#'       \item ps_nons: Estimated selection probabilities for non-probability sample
#'       \item variance_covariance: Variance-covariance matrix (if boot=FALSE)
#'       \item other diagnostic information
#'     }
#'   \item outcome: List of outcome model results including:
#'     \itemize{
#'       \item coefficients: Estimated outcome model coefficients
#'       \item std_err: Standard errors (if boot=FALSE)
#'       \item fitted.values: Fitted values for non-probability sample
#'       \item other model diagnostics
#'     }
#' }
#'
#' @details
#' The function implements a bias correction method for non-probability surveys
#' using both probability and non-probability samples. It fits both selection
#' and outcome models using maximum likelihood estimation via the nleqslv package.
#'
#' @note
#' The function may produce warnings if the optimization algorithm doesn't converge
#' or encounters numerical issues.
#'
#' @importFrom nleqslv nleqslv
#' @importFrom Matrix Diagonal
#'
#' @examples
#' \dontrun{
#' # Example code to be added
#' }
#'
#' @keywords models survey
#' @noRd
mm <- function(X,
               y,
               weights,
               weights_rand,
               R,
               n_nons,
               n_rand,
               method_selection,
               family,
               start_selection,
               start_outcome,
               maxit,
               nleqslv_method,
               nleqslv_global,
               nleqslv_xscalm,
               boot = FALSE) {

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  inv_link <- method$make_link_inv
  dinv_link <- method$make_link_inv_der

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  start <- c(start_outcome, start_selection) # TODO consider add info/error for end-user if one of starts provided only

  p <- ncol(X)
  if (is.null(start)) { # TODO add default start
    par0 <- rep(0, 2 * p)
  } else {
    par0 <- start
  }
  prior_weights <- c(weights_rand, weights)
  ######### NLESQLV
  multiroot <- nleqslv::nleqslv(
    x = par0,
    fn = u_theta_beta_dr,
    method = nleqslv_method,
    global = nleqslv_global,
    xscalm = nleqslv_xscalm,
    jacobian = TRUE,
    control = list(
      scalex = rep(1, length(par0)),
      maxit = maxit
    ), # TODO algorithm did not converge in maxit iterations for cloglog
    R = R,
    X = X,
    y = y,
    weights = prior_weights,
    method_selection = method_selection,
    family_nonprobsvy = family
  )
  par_sel <- multiroot$x
  if (multiroot$termcd %in% c(2:7, -10)) {
    switch(as.character(multiroot$termcd),
           "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv,
                         but user must check if function values are acceptably small."),
           "3" = warning("Algorithm did not find suitable point - has stalled cannot find
                         an acceptable new point when fitting selection model by nleqslv."),
           "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
           "5" = warning("Ill-conditioned Jacobian when fitting selection model by nleqslv."),
           "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
           "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
           "-10" = warning("User specified Jacobian is incorrect when fitting selection model by nleqslv.")
    )
  }

  theta_hat <- par_sel[1:(p)]
  beta_hat <- par_sel[(p + 1):(2 * p)]
  names(theta_hat) <- names(beta_hat) <- colnames(X)
  df_residual <- nrow(X) - length(theta_hat)

  # selection parameters
  ps <- inv_link(theta_hat %*% t(X)) # inv_link(as.vector(X_design %*% as.matrix(theta_hat)))
  eta_sel <- theta_hat %*% t(X)
  ps_der <- dinv_link(eta_sel)
  ps_nons <- ps[loc_nons]
  est_ps_rand <- ps[loc_rand]
  ps_nons_der <- ps_der[loc_nons]
  weights_nons <- 1 / ps_nons
  resids <- R - c(est_ps_rand, ps_nons)
  variance <- as.vector((t(resids) %*% resids) / df_residual)

  eta_out <- as.vector(beta_hat %*% t(X))
  y_hat <- family$linkinv(eta_out)
  y_rand_pred <- y_hat[loc_rand]
  y_nons_pred <- y_hat[loc_nons]

  if (!boot) {
    N_nons <- sum(weights * weights_nons)
    # variance-covariance matrix for selection model toFix
    V <- Matrix::Diagonal(n = length(ps), x = ps * (1 - ps))
    vcov_selection <- solve(t(X) %*% V %*% X)
    # vcov_selection <- matrix(0, nrow = nrow(X_design), ncol = ncol(X_design))
    theta_errors <- sqrt(diag(vcov_selection))
    # sigma_nons <- family$variance(mu = y_nons_pred, y = y[loc_nons])
    # sigma_rand <- family$variance(mu = y_rand_pred, y = y[loc_rand])
    sigma_nons <- family$variance(mu = y_nons_pred)
    sigma_rand <- family$variance(mu = y_rand_pred)
    residuals <- family$residuals(mu = y_nons_pred, y = y[loc_nons])
    # variance-covariance matrix for outcome model
    # vcov_outcome <- solve(t(X_design) %*% diag(sigma) %*% X_design)
    vcov_outcome <- solve(t(X[loc_nons, ]) %*% (sigma_nons * X[loc_nons, ]))
    beta_errors <- sqrt(diag(vcov_outcome))
  }




  # grad = multiroot$f.root[(p+1):(2*p)]
  if (!boot) {
    hess <- NA
    selection <- list(
      theta_hat = theta_hat, # TODO list as close as possible to SelecttionList
      grad = multiroot$fvec[1:(p)],
      hess = hess, # TODO
      ps_nons = ps_nons,
      est_ps_rand = est_ps_rand,
      variance_covariance = vcov_selection,
      df_residual = df_residual,
      log_likelihood = NA,
      eta = eta_sel,
      aic = NA,
      residuals = resids,
      variance = variance,
      method = method
    )

    outcome <- list(
      coefficients = beta_hat, # TODO list as close as possible to glm
      std_err = beta_errors,
      variance_covariance = vcov_outcome,
      df_residual = df_residual,
      family = list(
        mu = y_nons_pred,
        variance = sigma_nons,
        family = family$family
      ),
      residuals = residuals,
      fitted.values = y_nons_pred,
      sigma_rand = sigma_rand,
      y_rand_pred = y_rand_pred,
      y_nons_pred = y_nons_pred,
      linear.predictors = eta_out[loc_nons],
      X = X[loc_nons, ]
    )
  } else {
    selection <- list(
      coefficients = theta_hat, # TODO list as close as possible to SelectionList
      ps_nons = ps_nons
    )
    outcome <- list(
      coefficients = beta_hat,
      y_rand_pred = y_rand_pred, # TODO list as close as possible to SelectionList
      y_nons_pred = y_nons_pred
    )
  }

  list(
    selection = selection,
    outcome = outcome
  )
}
