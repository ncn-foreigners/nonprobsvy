#' @title Variance for Inverse Probability Weighted Estimator
#' @description Computes the variance for the inverse probability weighted (IPW) estimator
#' It is based on the selection method (`method_selection`) and estimation method `est_method`
#' @param svydesign Survey design object
#' @param X_nons Covariate matrix for non-sampled data
#' @param X_rand Covariate matrix for random sample
#' @param y_nons Response vector for non-sampled data
#' @param weights Weights for the observations
#' @param ps_nons Estimated propensity scores
#' @param mu_hat Estimated population mean
#' @param hess Hessian matrix
#' @param ps_nons_der Derivative of propensity scores for non-sampled data
#' @param N Population size
#' @param est_ps_rand Estimated propensity scores for random sample
#' @param ps_rand Probabilities of inclusion from the probability sample
#' @param est_ps_rand_der Derivative of estimated propensity scores for random sample
#' @param n_rand Number of observations in random sample
#' @param pop_size Population size
#' @param pop_totals Population totals
#' @param method_selection Method selection string
#' @param est_method Estimation method
#' @param theta Parameter vector
#' @param h type of the h() function
#' @param verbose Logical, if TRUE, prints additional information
#' @param var_cov1 Function to compute variance-covariance matrix for non-sampled data
#' @param var_cov2 Function to compute variance-covariance matrix for random sample
#' @return A list containing:
#' \describe{
#'   \item{var_nonprob}{Variance for non-probability component.}
#'   \item{var_prob}{Variance for probability component.}
#'   \item{var}{Total variance.}
#' }
#' @keywords internal
#' @noRd
internal_varIPW <- function(svydesign,
                            X_nons,
                            X_rand,
                            y_nons,
                            weights,
                            ps_nons,
                            mu_hat,
                            hess,
                            ps_nons_der,
                            N,
                            est_ps_rand,
                            ps_rand,
                            est_ps_rand_der,
                            n_rand,
                            pop_size,
                            pop_totals,
                            method_selection,
                            est_method,
                            theta,
                            gee_h_fun,
                            verbose,
                            var_cov1,
                            var_cov2) {

  eta <- as.vector(X_nons %*% as.matrix(theta))
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
  b_obj <- method$b_vec_ipw(
    X = X_nons,
    ps = ps_nons,
    psd = ps_nons_der,
    y = y_nons,
    mu = mu_hat,
    hess = hess,
    eta = eta,
    pop_size = pop_size,
    weights = weights,
    verbose = verbose
  )
  b <- b_obj$b

  # sparse matrix
  b_vec <- cbind(-1, b)
  H_mx <- try(cbind(0, N * solve(hess)), silent = TRUE)
  if (inherits(H_mx, "try-error")) {
    if (verbose) message("solve() failed, using ginv() instead.")
    H_mx <- cbind(0, N * ginv(hess))
  }
  sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

  V1 <- var_cov1(
    X = X_nons,
    y = y_nons,
    mu = mu_hat,
    ps = ps_nons,
    psd = ps_nons_der,
    pop_size = pop_size,
    est_method = est_method,
    gee_h_fun = gee_h_fun,
    weights = weights,
    pop_totals = pop_totals
  ) # fixed
  V2 <- var_cov2(
    X = X_rand,
    svydesign = svydesign,
    eps = est_ps_rand,
    est_method = est_method,
    gee_h_fun = gee_h_fun,
    pop_totals = pop_totals,
    psd = est_ps_rand_der
  )

  # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
  V_mx_nonprob <- sparse_mx %*% V1 %*% t(as.matrix(sparse_mx)) # nonprobability component
  V_mx_prob <- sparse_mx %*% V2 %*% t(as.matrix(sparse_mx)) # probability component
  V_mx <- V_mx_nonprob + V_mx_prob

  var_nonprob <- as.vector(V_mx_nonprob[1, 1])
  var_prob <- as.vector(V_mx_prob[1, 1])
  var <- as.vector(V_mx[1, 1])

  list(
    var_nonprob = var_nonprob,
    var_prob = var_prob,
    var = var
  )
}
