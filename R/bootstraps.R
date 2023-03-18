# These functions are only used internally, so there is no need for documenting them
#' @importFrom survey as.svrepdesign

bootMI <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   svydesign,
                   rep_type,
                   ...
                   ){


  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  N <- sum(weights_rand)
  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
  k <- 1

  while (k <= num_boot) {

    strap <- sample.int(replace = TRUE, n = n_nons)
    weights_strap <- weights[strap]
    X_nons_strap <- X_nons[strap,]
    y_strap <- y[strap]

    #using svy package
    strap_rand_svy <- which(rep_weights[,k] != 0)
    weights_rand_strap_svy <- rep_weights[,k] * weights_rand


    model_strap <- nonprobMI_fit(x = X_nons_strap,
                                 y = y_strap,
                                 weights = weights_strap,
                                 family_outcome = family_outcome)

    beta <- model_strap$coefficients

    ystrap_rand <- as.numeric(X_rand %*% beta)

    mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  boot_var

}

bootIPW <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    family_outcome,
                    num_boot,
                    weights_rand,
                    mu_hat,
                    dependency,
                    ...){





}

bootDR <- function(SelectionModel,
                   OutcomeModel,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   method_selection,
                   ...){

  method <- get_method(method_selection)

  ps_method <- method$make_propen_score # function for propensity score estimation
  loglike <- method$make_log_like
  gradient <- method$make_gradient
  hessian <- method$make_hessian

  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  N <- sum(weights_rand)
  k <- 1

  while (k <= num_boot) {

    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)
    weights_strap <- weights[strap_nons]
    X_nons_strap <- X_nons[strap_nons, ]
    y_strap <- y[strap_nons]

    X_rand <- X_rand[strap_rand, ]
    weights_rand_strap <- weights_rand[strap_rand]


    model_strap <- nonprobMI_fit(x = X_nons_strap,
                                 y = y_strap,
                                 weights = weights_strap,
                                 family_outcome = family_outcome)
    beta <- model_strap$coefficients
    ystrap_rand <- as.numeric(X_rand %*% beta)
    ystrap_nons <- as.numeric(X_nons %*% beta)

    log_like <- loglike(X_nons, X_rand, weights_rand)
    gradient <- gradient(X_nons, X_rand, weights_rand)
    hessian <- hessian(X_nons, X_rand, weights_rand)

    maxLik_nons_obj <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)
    maxLik_rand_obj <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)

    ps_nons <- maxLik_nons_obj$ps
    est_ps_rand <- maxLik_rand_obj$ps
    hess <- maxLik_nons_obj$hess
    theta_hat <- maxLik_nons_obj$theta_hat

    N_est_nons <- sum(1/ps_nons)
    N_est_rand <- sum(weights_rand_strap)

    mu_hat_boot <- mu_hatDR(y = y_strap,
                            y_nons = y_nons_pred,
                            y_rand = ystrap_rand,
                            weights_nons = weights_strap,
                            weights_rand = weights_rand_strap,
                            N_nons = N_est_nons,
                            N_rand = N_est_rand)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mean(mu_hats))^2)
  boot_var


}
