#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom sampling calib
#' @export
#' @rdname main_doc

nonprobCL <- function(outcome,
                      data,
                      svydesign,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_outcome = controlOut(),
                      control_inference = controlInf(var_method = "analytic"),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  y_name <- colnames(XY_nons)[1]
  X_nons <- model.matrix(XY_nons, data) # matrix of the  nonprobability sample
  svydesign$variables[,y_name] <- rep(0, nrow(svydesign$variables))
  X_rand <- model.matrix(outcome, svydesign$variables) # matrix of the probability sample
  y_nons <- XY_nons[,1]

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand)

  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand
  N_est_rand <- sum(weights_rand)

  total <- t(X_rand) %*% as.matrix(weights_rand)
  weights_nons <- sampling::calib(Xs = X_nons, d = weights, total = total, method = "linear")
  N_est_nons <- sum(weights_nons)
  mu_hat <- 1/N_est_nons * sum(weights_nons * y_nons)

  structure(
    list(mean = mu_hat),
    class = "Calibration Weighting")
}
