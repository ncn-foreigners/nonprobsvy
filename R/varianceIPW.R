# Variance for inverse probability weighted estimator
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
                            h,
                            verbose,
                            var_cov1 = var_cov1,
                            var_cov2 = var_cov2) {
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
    h = h,
    weights = weights,
    pop_totals = pop_totals
  ) # fixed
  V2 <- var_cov2(
    X = X_rand,
    svydesign = svydesign,
    eps = est_ps_rand,
    est_method = est_method,
    h = h,
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
  # vector of variances for theta_hat
  # theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

  list(
    var_nonprob = var_nonprob,
    var_prob = var_prob,
    var = var
  )
}
