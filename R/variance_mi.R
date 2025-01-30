# Variance for mass imputation estimator
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom stats loess.control
#' @importFrom survey svymean
internal_varMI <- function(svydesign,
                           X_nons,
                           X_rand,
                           y,
                           y_pred,
                           y_hat,
                           weights_rand,
                           method,
                           n_rand,
                           n_nons,
                           N,
                           family,
                           model_obj,
                           pop_totals,
                           k,
                           predictive_match,
                           nn_exact_se,
                           pmm_reg_engine,
                           pi_ij) {
  parameters <- model_obj$parameters

  if (is.character(family)) {
    family_nonprobsvy <- paste(family, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (is.null(pop_totals)) {
    svydesign_mean <- survey::svymean(~y_hat_MI, svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var")) # probability component, should be bigger for nn
    if (method == "nn") {
      sigma_hat <- mean((y - y_pred)^2) # family_nonprobsvy$variance(mu = y_pred, y  = y)
      est_ps <- n_nons / N
      var_nonprob <- n_rand / N^2 * (1 - est_ps) / est_ps * sigma_hat

      if (nn_exact_se) {
        var_nonprob <- nn_exact(
          pi_ij        = pi_ij,
          weights_rand = weights_rand,
          n_nons       = n_nons,
          y            = y,
          X_nons       = X_nons,
          X_rand       = X_rand,
          k            = k,
          # TODO:: add control here
          # control      = control
          N            = N
        )
      }
    } else if (method == "glm") { # TODO add variance for count binary outcome variable control_outcome$method

      beta <- parameters[, 1]
      eta_nons <- X_nons %*% beta
      eta_rand <- X_rand %*% beta

      mx <- 1 / N * colSums(as.data.frame(X_rand) * (weights_rand * family_nonprobsvy$mu.eta(eta_rand)))
      c <- solve(1 / n_nons * t(as.data.frame(X_nons) * family_nonprobsvy$mu.eta(eta_nons)) %*% X_nons) %*% mx
      residuals <- family_nonprobsvy$residuals(mu = y_pred, y = y)

      # nonprobability component
      var_nonprob <- 1 / n_nons^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2
      var_nonprob <- as.vector(var_nonprob)
    } else if (method == "pmm") {
      var_prob <- as.numeric(attr(svymean(~y_hat_MI, svydesign), "var"))

      # This in general cannot be computed from sample itself, we need to make
      # a bootstrap. Sometimes this term is negligible hence by default its
      # not computed, but it should be computed in serious publications
      var_nonprob <- 0

      # An option in controlInf controls this
      # Maybe add a warning/message if this computation is omited
      if (nn_exact_se) {
        var_nonprob <- pmm_exact(
          pi_ij = pi_ij,
          weights_rand = weights_rand,
          n_nons = n_nons,
          y = y,
          pmm_reg_engine = pmm_reg_engine,
          model_obj = model_obj,
          svydesign = svydesign,
          predictive_match = predictive_match,
          k = k,
          N = N
        )
      }
    }
  } else {
    if (method == "nn") {
      sigma_hat <- mean((y - y_pred)^2) # family_nonprobsvy$variance(mu = y_pred, y  = y)
      est_ps <- n_nons / N
      var_nonprob <- n_nons / N^2 * (1 - est_ps) / est_ps * sigma_hat # what instead of n_rand here (?) now just n_nons
    } else if (method == "glm") {
      beta <- parameters[, 1]
      eta_nons <- X_nons %*% beta
      if (family %in% c("binomial", "poisson")) { # TODO consider this chunk of code
        eta_rand <- pop_totals %*% beta / pop_totals[1]
      } else {
        eta_rand <- pop_totals %*% beta
      }
      mx <- 1 / N * pop_totals * as.vector(family_nonprobsvy$mu.eta(eta_rand))
      c <- solve(1 / n_nons * t(as.data.frame(X_nons) * family_nonprobsvy$mu.eta(eta_nons)) %*% X_nons) %*% mx
      residuals <- family_nonprobsvy$residuals(mu = y_pred, y = y)

      # nonprobability component
      var_nonprob <- 1 / n_nons^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2
      var_nonprob <- as.vector(var_nonprob)
    } else if (method == "pmm") {
      # TODO
      # beta <- parameters[,1]
      # eta_nons <- X_nons %*% beta
      #
      # if (family %in% c("binomial", "poisson")) { # TODO consider this chunk of code
      #   eta_rand <- pop_totals %*% beta / pop_totals[1]
      # } else {
      #   eta_rand <- pop_totals %*% beta
      # }
      #
      # residuals <- family_nonprobsvy$residuals(mu = y_pred, y  = y)

      # nonprobability component
      # var_nonprob <- 1/n_nons^2 * t(as.matrix(residuals^2)) %*% (family_nonprobsvy$mu_der(eta_nons) %*% t(X_nons))^2
      var_nonprob <- 0
      var_nonprob <- as.vector(var_nonprob)
    }
    var_prob <- 0
  }
  list(
    var_prob = var_prob,
    var_nonprob = var_nonprob
  )
}
