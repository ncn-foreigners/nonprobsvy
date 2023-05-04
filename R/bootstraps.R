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
                    R,
                    family_outcome,
                    num_boot,
                    weights_rand,
                    mu_hat,
                    method_selection,
                    n_nons,
                    n_rand,
                    optim_method,
                    est_method,
                    h,
                    maxit,
                    pop_size = NULL,
                    ...){
  mu_hats <- vector(mode = "numeric", length = num_boot)
  N <- sum(weights_rand)
  k <- 1
  est_method_fun <- get(est_method, mode = "function", envir = parent.frame())

  while (k <= num_boot) {

    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)

    X <- rbind(X_rand[strap_rand, ],
               X_nons[strap_nons, ])

    model_sel <- internal_selection(X = X,
                                    X_nons = X_nons[strap_nons, ],
                                    X_rand = X_rand[strap_rand, ],
                                    weights = weights[strap_nons],
                                    weights_rand = weights_rand[strap_rand],
                                    R = R,
                                    method_selection = method_selection,
                                    optim_method = optim_method,
                                    h = h,
                                    est_method =  est_method,
                                    maxit = maxit)
    est_method_obj <- est_method_fun(model = model_sel,
                                     method_selection = method_selection)

    ps_nons <- est_method_obj$ps_nons
    weights_nons <- 1/ps_nons
    N_est_nons <- ifelse(is.null(pop_size), sum(1/ps_nons), pop_size)

    mu_hat_boot <- mu_hatIPW(y = y[strap_nons],
                            weights = weights_nons,
                            N = N_est_nons) # IPW estimator
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  boot_var
}

bootDR <- function(SelectionModel,
                   OutcomeModel,
                   family_outcome,
                   num_boot,
                   weights,
                   weights_rand,
                   R,
                   mu_hat,
                   method_selection,
                   n_nons,
                   n_rand,
                   optim_method,
                   est_method,
                   h,
                   maxit = NULL,
                   ...) {

  mu_hats <- vector(mode = "numeric", length = num_boot)
  N <- sum(weights_rand)
  k <- 1
  est_method_fun <- get(est_method, mode = "function", envir = parent.frame())

  while (k <= num_boot) {

    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)

    model_out <- internal_outcome(X_nons = OutcomeModel$X_nons[strap_nons, ],
                                  X_rand = OutcomeModel$X_rand[strap_rand, ],
                                  y = OutcomeModel$y[strap_nons],
                                  weights = weights[strap_nons],
                                  family_outcome = family_outcome)

    y_rand_pred <- model_out$y_rand_pred
    y_nons_pred <- model_out$y_nons_pred
    model_nons_coefs <- model_out$model_nons_coefs

    X_sel <- rbind(SelectionModel$X_rand[strap_rand, ],
                   SelectionModel$X_nons[strap_nons, ])

    model_sel <- internal_selection(X = X_sel,
                                    X_nons = SelectionModel$X_nons[strap_nons, ],
                                    X_rand = SelectionModel$X_rand[strap_rand, ],
                                    weights = weights[strap_nons],
                                    weights_rand = weights_rand[strap_rand],
                                    R = R,
                                    method_selection = method_selection,
                                    optim_method = optim_method,
                                    h = h,
                                    est_method = est_method,
                                    maxit = maxit)

    est_method_obj <- est_method_fun(model = model_sel,
                                      method_selection = method_selection)
    ps_nons <- est_method_obj$ps_nons
    weights_nons <- 1/ps_nons
    N_est_nons <- sum(weights_nons)
    N_est_rand <- sum(weights_rand[strap_rand])

    mu_hat_boot <- mu_hatDR(y = OutcomeModel$y_nons[strap_nons],
                            y_nons = y_nons_pred,
                            y_rand = y_rand_pred,
                            weights_nons = weights_nons,
                            weights_rand = weights_rand[strap_rand],
                            N_nons = N_est_nons,
                            N_rand = N_est_rand)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  boot_var
}
