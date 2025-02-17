# Variance for doubly robust estimator
internal_varDR <- function(outcome_model,
                           selection_model,
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
                           gee_h_fun,
                           pop_totals,
                           sigma,
                           bias_correction,
                           verbose) {
  ######### mm
  if (bias_correction == TRUE) {
    infl1 <- (weights * (outcome_model$y_nons - y_nons_pred))^2 / ps_nons^2
    infl2 <- (weights * (outcome_model$y_nons - y_nons_pred))^2 / ps_nons

    # Variance estimators ####
    svydesign <- stats::update(svydesign,
      y_rand = y_rand_pred
    )
    svydesign_mean <- survey::svymean(~y_rand, svydesign)

    var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
    var_nonprob <- (sum((infl1) - 2 * infl2) + sum(weights_rand * sigma)) / N_nons^2 # TODO potential bug here nonprobability component
  } else {
    eta <- as.vector(selection_model$X_nons %*% as.matrix(theta))
    h_n <- 1 / N_nons * sum(outcome_model$y_nons - y_nons_pred) # TODO add weights # errors mean

    method <- switch(method_selection,
                     "logit" = model_ps("logit"),
                     "probit" = model_ps("probit"),
                     "cloglog" = model_ps("cloglog"))

    est_method <- get_method(est_method)

    b <- method$b_vec_dr(
      X = selection_model$X_nons,
      ps = ps_nons,
      psd = ps_nons_der,
      y = outcome_model$y_nons,
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
      y = outcome_model$y_nons,
      y_pred = y_nons_pred,
      h_n = h_n,
      X = selection_model$X_nons,
      b = b,
      N = N_nons,
      gee_h_fun = gee_h_fun,
      method_selection = method_selection,
      weights = weights,
      pop_totals = pop_totals
    )

    if (is.null(pop_totals)) {
      t_comp <- est_method$make_t_comp(
        X = selection_model$X_rand,
        ps = est_ps_rand,
        psd = est_ps_rand_der,
        b = b,
        gee_h_fun = gee_h_fun,
        y_rand = y_rand_pred,
        y_nons = y_nons_pred,
        N = N_nons,
        method_selection = method_selection,
        weights = weights
      )
      # design based variance estimation based on approximations of the second-order inclusion probabilities
      svydesign <- stats::update(svydesign, t_comp = t_comp)
      svydesign_mean <- survey::svymean(~t_comp, svydesign) # perhaps using survey package to compute prob variance
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
