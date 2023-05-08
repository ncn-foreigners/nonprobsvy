
# Object with output parameters for Maximum likelihood estimation for propensity scores
mle <- function(...) {

  estimation_model <- function(model, method_selection, ...) {
    maxLik_nons_obj <- model$maxLik_nons_obj
    maxLik_rand_obj <- model$maxLik_rand_obj
    log_likelihood <- model$log_likelihood # maximum of the loglikelihood function
    theta_hat <- model$theta

    ps_nons <- maxLik_nons_obj$ps
    est_ps_rand <- maxLik_rand_obj$ps
    hess <- maxLik_nons_obj$hess
    grad <- maxLik_rand_obj$grad
    var_cov1 <- model$var_cov1
    var_cov2 <- model$var_cov2
    ps_nons_der <- NULL
    est_ps_rand_der <- NULL

    if (method_selection == "probit") { # for probit model, propensity score derivative is required
      ps_nons_der <- maxLik_nons_obj$psd
      est_ps_rand_der <- maxLik_rand_obj$psd
    }

    list(theta_hat = theta_hat,
         grad = grad,
         hess = hess,
         var_cov1 = var_cov1,
         var_cov2 = var_cov2,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ps_nons_der,
         est_ps_rand_der = est_ps_rand_der)
  }

  make_t <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection) {
    method <- get_method(method_selection)
    t <- method$t_vec(X = SelectionModel$X_rand,
                      ps = est_ps_rand,
                      psd = est_ps_rand_der,
                      b = b,
                      y_rand = y_rand_pred,
                      y_nons = y_nons_pred,
                      N = N_nons)
    t
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection) {
    method <- get_method(method_selection)
    var_nonprob <-  method$var_nonprob(ps = ps_nons,
                                       psd = ps_nons_der,
                                       y = OutcomeModel$y_nons,
                                       y_pred = y_nons_pred,
                                       h_n = h_n,
                                       X = SelectionModel$X_nons,
                                       b = b,
                                       N = N_nons)
    as.numeric(var_nonprob)

  }

  structure(
    list(estimation_model = estimation_model,
         make_t = make_t,
         make_var_nonprob = make_var_nonprob),
    class = "est_method"
  )

}
# Object with output parameters for estimation by Generalized Estimating Equations for propensity scores
gee <- function(...) {

  estimation_model <- function(model, method_selection) {
    theta_hat <- model$theta_hat
    hess <- model$hess
    grad <- model$grad
    ps_nons <- model$ps_nons
    est_ps_rand <- model$est_ps_rand
    ps_nons_der <- model$ps_nons_der
    est_ps_rand_der <- model$est_ps_rand_der
    var_cov1 <- model$var_cov1
    var_cov2 <- model$var_cov2

    list(theta_hat = theta_hat,
         grad = grad,
         hess = hess,
         var_cov1 = var_cov1,
         var_cov2 = var_cov2,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ps_nons_der,
         est_ps_rand_der = est_ps_rand_der)
  }

  make_t <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection) {
    if (h == "1") {
      t <- X %*% t(as.matrix(b)) + y_rand - 1/N * sum(y_nons)
    } else if (h == "2") {
      t <- as.vector(ps) * X %*% t(as.matrix(b)) + y_rand - 1/N * sum(y_nons)
    }
    t
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection) {
    if (h == "1") {
      var_nonprob <- 1/N^2 * sum((1 - ps) * (((OutcomeModel$y - y_pred - h_n)/ps_nons) - b %*% t(X))^2)
    } else if (h == "2") {
      var_nonprob <- 1/N^2 * sum((1 - ps) * (((OutcomeModel$y - y_pred - h_n) - b %*% t(X))/ps)^2)
    }
    as.numeric(var_nonprob)
  }

  structure(
    list(estimation_model = estimation_model,
         make_t = make_t,
         make_var_nonprob = make_var_nonprob),
    class = "est_method"
  )

}
