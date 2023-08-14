# Internal functions, no need for documenting them
nonprobOv <- function(X_nons,
                      X_rand,
                      weights,
                      weights_rand,
                      dependent,
                      method_selection,
                      key_var_prob,
                      control,
                      idx_nonprob,
                      idx_prob,
                      ov_method) {

  betamodel <- betareg::betareg.fit(x = X_rand,
                                    y = 1/weights_rand,
                                    link = method_selection) # Assumption A3 from the paper

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv
  prob <- inv_link(X_nons %*% betamodel$coefficients$mean)
  weights_rnons <- 1/prob

  # weights_rnons are inclusion probabilities in the probability sample for units in the nonprobability sample
  weights_rr <- c(weights_rnons, weights_rand)

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X <- rbind(X_nons, X_rand)
  prior_weights <- c(weights[-idx_nonprob], weights_rand[-idx_prob])
  XX <- rbind(X_nons[-idx_nonprob,], X_rand[-idx_prob,])
  RR_nons <- rep(1, nrow(X_nons[-idx_nonprob,]))
  RR_rand <- rep(0, nrow(X_rand[-idx_prob,]))
  RR <- c(RR_nons, RR_rand)
  locc_nons <- which(RR == 1)
  locc_rand <- which(RR == 0)

  #### MLE
  if (control$ov_method == "mle") {
    start <- start_fit(X = XX,
                       R = RR,
                       weights = prior_weights[locc_nons],
                       weights_rand = prior_weights[locc_rand],
                       method_selection = method_selection)
    model <- method$make_max_lik(X_nons = XX[locc_nons,],
                                 X_rand = XX[locc_rand,],
                                 prior_weights[locc_nons],
                                 prior_weights[locc_rand],
                                 start = start,
                                 control = control)
    eta <- X %*% model$theta_hat
    ps <- inv_link(eta)
    var_cov <- solve(-model$hess)
    theta_errors <- sqrt(diag(var_cov))
    parameters <- matrix(c(model$theta_hat, theta_errors), # TODO
                         ncol = 2,
                         dimnames = list(names(XX),
                                         c("Estimate", "Std. Error")))
    L <- NULL
    O <- NULL
    # TODO add summary

    if (dependent) {
      stop("TODO")
    }
  } else if (control$ov_method == "gee") {
    ### GEE
    # TODO problem with hess/jac
    model <- theta_h_estimation(R = RR,
                                X = XX,
                                weights_rand = prior_weights[locc_rand],
                                weights = prior_weights[locc_nons],
                                h = control$h_x,
                                method_selection,
                                maxit = control$maxit,
                                pop_totals = NULL,
                                pop_means = NULL)
    eta <- X %*% model$theta_h
    var_cov <- solve(-model$hess)
    theta_errors <- sqrt(diag(var_cov))
    ps <- inv_link(eta)
    parameters <- matrix(c(model$theta_h, theta_errors), # TODO
                         ncol = 2,
                         dimnames = list(names(XX),
                                         c("Estimate", "Std. Error")))
    L <- NULL
    O <- NULL
    if (dependent) {
      stop("TODO")
    }
  } else if (control$ov_method == "ev") {
    O_hat <- O_hat_model(x = XX,
                         X = X,
                         R = RR,
                         prior_weights = prior_weights,
                         method_selection = method_selection,
                         method_obj = method,
                         control = control,
                         loc_nons = locc_nons,
                         loc_rand = locc_rand)
    O <- O_hat$O
    parameters <- O_hat$parameters
    L <- NULL
    if (dependent) {
      L <- L_hat_model(x = X_rand,
                       X = X,
                       key_var = key_var_prob,
                       method_selection = method_selection) # key instead of R_rand here

      ps <- O * L/(weights_rr - 1)
    } else {
      ps <- O/(O + weights_rr - 1)
    }
  } else {
    stop("Invalid method for overlap model.")
  }

  structure(
    list(est_ps_rand <- ps[loc_rand],
         ps_nons = ps[loc_nons],
         O_hat = O[loc_nons],
         L_hat = L,
         weights_rnons = weights_rnons,
         parameters = parameters)
  )
}

boot_overlap <- function(X_rand,
                         X_nons,
                         y,
                         weights_nons, # for the nonprobability sample
                         weights,
                         mu_hat,
                         O_hat,
                         L_hat,
                         weights_rand,
                         weights_rnons,
                         method_selection,
                         dependency,
                         N,
                         type,
                         idx_nonprob,
                         idx_prob,
                         control,
                         family_outcome = "gaussian",
                         num_boot = 1000) {


  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  mu_hats <- vector(mode = "numeric", length = num_boot)
  weights_rr <- c(weights_rand, weights_rnons)
  R_nons <- rep(1, n_nons)
  R_rand <- rep(0, n_rand)
  R <- c(R_nons, R_rand)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }


  weights_nons <- weights_nons*N/sum(weights_nons)
  prob_floor <- weights_nons - floor(weights_nons)
  weights_nons_strap <- weights_nons
  weights_idx <- sample.int(length(weights_nons), replace = TRUE, prob = prob_floor)
  weights_nons_strap[weights_idx] <- ceiling(weights_nons_strap[weights_idx])
  weights_nons_strap[!weights_idx] <- floor(weights_nons_strap[!weights_idx])

  pseudo_pop <- X_nons[rep(1:n_nons, times = weights_nons_strap), ] # pseudo population - to consider
  pseudo_y <- y[rep(1:n_nons, times = weights_nons_strap)]
  pseudo_weights_rnons <- rep(weights_rnons, times = weights_nons_strap)
  pseudo_weights <- rep(weights_nons, times = weights_nons_strap)


  k <- 1

  if (!dependency) {
    while (k <= num_boot) {

      strap_rand <- sample.int(n = nrow(pseudo_pop), size = n_rand, replace = TRUE, prob = 1/pseudo_weights_rnons)
      strap_nons <- sample.int(n = nrow(pseudo_pop), size = n_nons, replace = TRUE, prob = 1/pseudo_weights)
      #strap_d <- sample.int(n = n_nons + n_rand, replace = TRUE)
      weights_rnons_strap <- pseudo_weights_rnons[c(strap_rand, strap_nons)]

      X_rand_strap <- pseudo_pop[strap_rand, ]
      X_nons_strap <- pseudo_pop[strap_nons, ]
      X_strap <- rbind(X_nons_strap, X_rand_strap)
      XX_strap <- rbind(X_nons_strap[-idx_nonprob,], X_rand_strap[-idx_prob,])
      R_nons_strap <- rep(1, nrow(X_nons_strap))
      R_rand_strap <- rep(0, nrow(X_rand_strap))
      R_strap <- c(R_nons_strap, R_rand_strap)
      RR_nons_strap <- rep(1, nrow(X_nons_strap[-idx_nonprob,]))
      RR_rand_strap <- rep(0, nrow(X_rand_strap[-idx_prob,]))
      RR_strap <- c(RR_nons_strap, RR_rand_strap)
      loc_nons_strap <- which(R_strap == 1)
      loc_rand_strap <- which(R_strap == 0)
      locc_nons_strap <- which(RR_strap == 1)
      locc_rand_strap <- which(RR_strap == 0)

      O_strap_model <- O_hat_model(x = XX_strap,
                                   X = X_strap,
                                   R = RR_strap,
                                   prior_weights = c(weights[-idx_nonprob], weights_rand[-idx_prob]), # TODO to consider
                                   method_selection = method_selection,
                                   method_obj = method,
                                   control = control,
                                   loc_nons = locc_nons_strap,
                                   loc_rand = locc_rand_strap
                                  )
      O_strap <- O_strap_model$O

      ps_strap <- O_strap/(O_strap + weights_rnons_strap - 1)
      weights_strap <- 1/ps_strap[loc_nons_strap]
      N_strap <- sum(weights_strap)

      if (type == "IPW") {
      mu_hat_boot <- mu_hatIPW(y = pseudo_y[strap_nons],
                               weights = weights, # TODO
                               weights_nons = weights_strap,
                               N = N_strap)
      } else if (type == "DR") {
        if(is.character(family_outcome)) {
          family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
          family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
          family_nonprobsvy <- family_nonprobsvy()
        }
        model <- stats::glm.fit(x = X_nons_strap,
                                y = pseudo_y[strap_nons],
                                weights = weights,
                                family = get_method(family_outcome),
                                # control = list(control$epsilon,
                                #                control$maxit,
                                #                control$trace),
                                intercept = FALSE)
        model_nons_coefs <- model$coefficients
        eta <- X_rand_strap %*% model_nons_coefs
        y_rand_pred <- family_nonprobsvy$mu(eta)
        y_nons_pred <- model$fitted.values

        N_rand <- sum(weights_rand)
        mu_hat_boot <- mu_hatDR(y = pseudo_y[strap_nons],
                                y_nons = y_nons_pred,
                                y_rand = y_rand_pred,
                                weights = weights, # TODO
                                weights_nons = weights_strap,
                                weights_rand = weights_rand, #TODO
                                N_nons = N_strap,
                                N_rand = N_rand)
      } else {
        stop("Invalid population mean estimator type")
      }
      mu_hats[k] <- mu_hat_boot
      k <- k + 1
    }

  } else {

    pseudo_L <- rep(L_hat[loc_nons], times = weights_nons_strap)
    pseudo_O <- rep(O_hat[loc_nons], times = weights_nons_strap)

    ps1 <- 1 - pseudo_L
    ps0 <- 1 - pseudo_weights_rnons + pseudo_weights_rnons*pseudo_O*pseudo_L + pseudo_L*(pseudo_weights_rnons - 1)


    # weights for strap_nons to compute

    while (k <= num_boot) {

      strap_rand <- sample.int(n = nrow(pseudo_pop), size = n_rand, replace = TRUE, prob = 1/pseudo_weights_rnons)

      ps1[!strap_rand] <- ps0
      strap_nons <- sample.int(n = nrow(pseudo_pop), size = n_nons, replace = TRUE, prob = 1/ps1)
      weights_rnons_strap <- pseudo_weights_rnons[c(strap_rand, strap_nons)]

      X_rand_strap <- pseudo_pop[strap_rand, ]
      X_nons_strap <- pseudo_pop[strap_nons, ]
      X_strap <- rbind(X_nons_strap, X_rand_strap)
      R_nons_strap <- rep(1, nrow(X_nons_strap))
      R_rand_strap <- rep(0, nrow(X_rand_strap))
      R_strap <- c(R_nons_strap, R_rand_strap)
      loc_nons_strap <- which(R_strap == 1)
      loc_rand_strap <- which(R_strap == 0)
      # TODO Removing overlapping units here (?)

      O_strap_model <- O_hat_model(x = X_strap,
                                   X = X_strap,
                                   R = R_strap,
                                   prior_weights = c(weights, weights_rand), # TODO to consider
                                   method_selection = method_selection,
                                   control = control,
                                   loc_nons = loc_nons_strap,
                                   loc_rand = loc_rand_strap)
      O_strap <- O_strap_model$O

      L_strap <- L_hat_model(x = X_rand_strap,
                             X = X_strap,
                             key_var = R_rand_strap, # TODO
                             method_selection = method_selection)

      ps_strap <- O_strap * L_strap/(weights_rnons_strap - 1)
      weights_strap <- 1/ps_strap[loc_nons_strap]
      N_strap <- sum(weights_strap)

      mu_hat_boot <- mu_hatIPW(y = pseudo_y[strap_nons],
                               weights = weights,
                               weights_nons = weights_nons_strap,
                               N = N_strap)
      mu_hats[k] <- mu_hat_boot
      k <- k + 1
    }
  }
  var <- 1/(num_boot-1) * sum((mu_hats - mu_hat)^2)
  var
}


L_hat_model <- function(x,
                        X,
                        key_var,
                        method_selection) {
  modelL <- stats::glm.fit(x = x,
                           y = key_var,
                           family = binomial(link = method_selection))
  L <- 1/(1 + exp(X %*% modelL$coefficients))
  L
}


O_hat_model <- function(x,
                        X,
                        R,
                        prior_weights,
                        method_selection,
                        method_obj,
                        control,
                        loc_nons,
                        loc_rand) { # add method_selection
    modelO <- stats::glm.fit(x = x,
                             y = R,
                             family = binomial(link = method_selection),
                             # TODO weights = prior_weights
                             )
    O <- exp(X %*% modelO$coefficients)
    parameters = summary.glm(modelO)$coefficients

  list(O = O,
       parameters = parameters)
}
