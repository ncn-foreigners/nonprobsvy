# Variance for doubly robust estimator
# TODO add nn and pmm
internal_varDR <- function(OutcomeModel,
                           SelectionModel,
                           y_nons_pred,
                           weights,
                           weights_rand,
                           method_selection,
                           control_selection,
                           theta,
                           ps_nons,
                           hess,
                           ps_nons_der,
                           est_ps_rand,
                           y_rand_pred,
                           N_nons,
                           est_ps_rand_der,
                           svydesign,
                           est_method,
                           h,
                           pop_totals,
                           sigma,
                           bias_correction,
                           verbose) {
  ######### mm
  if (bias_correction == TRUE) {
    infl1 <- (weights * (OutcomeModel$y_nons - y_nons_pred))^2 / ps_nons^2
    infl2 <- (weights * (OutcomeModel$y_nons - y_nons_pred))^2 / ps_nons

    # Variance estimators ####
    svydesign <- stats::update(svydesign,
      y_rand = y_rand_pred
    )
    svydesign_mean <- survey::svymean(~y_rand, svydesign)

    var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
    var_nonprob <- (sum((infl1) - 2 * infl2) + sum(weights_rand * sigma)) / N_nons^2 # TODO potential bug here nonprobability component
  } else {
    eta <- as.vector(SelectionModel$X_nons %*% as.matrix(theta))
    h_n <- 1 / N_nons * sum(OutcomeModel$y_nons - y_nons_pred) # TODO add weights # errors mean
    method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
    method <- get_method(method_selection)
    est_method <- get_method(est_method)

    b <- method$b_vec_dr(
      X = SelectionModel$X_nons,
      ps = ps_nons,
      psd = ps_nons_der,
      y = OutcomeModel$y_nons,
      hess = hess,
      eta = eta,
      h_n = h_n,
      y_pred = y_nons_pred,
      weights = weights,
      verbose = verbose
    )

    # asymptotic variance by each propensity score method (nonprobability component)
    var_nonprob <- est_method$make_var_nonprob(
      ps = ps_nons,
      psd = ps_nons_der,
      y = OutcomeModel$y_nons,
      y_pred = y_nons_pred,
      h_n = h_n,
      X = SelectionModel$X_nons,
      b = b,
      N = N_nons,
      h = h,
      method_selection = method_selection,
      weights = weights,
      pop_totals = pop_totals
    )

    if (is.null(pop_totals)) {
      t <- est_method$make_t(
        X = SelectionModel$X_rand,
        ps = est_ps_rand,
        psd = est_ps_rand_der,
        b = b,
        h = h,
        y_rand = y_rand_pred,
        y_nons = y_nons_pred,
        N = N_nons,
        method_selection = method_selection,
        weights = weights
      )
      # design based variance estimation based on approximations of the second-order inclusion probabilities
      svydesign <- stats::update(svydesign,
        t = t
      )
      svydesign_mean <- survey::svymean(~t, svydesign) # perhaps using survey package to compute prob variance
      var_prob <- as.vector(attr(svydesign_mean, "var"))
    } else {
      var_prob <- 0
    }
  }
  list(
    var_prob = var_prob,
    var_nonprob = var_nonprob
  )
}
