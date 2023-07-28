# These functions are only used internally, so there is no need for documenting them
#' @importFrom survey as.svrepdesign
#' @importFrom nleqslv nleqslv

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
                   method,
                   control,
                   ...
                   ){ # TODO add methods instead of conditional loops

  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  N <- sum(weights_rand)
  k <- 1
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (method == "glm") {
    rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

    while (k <= num_boot) {

      strap <- sample.int(replace = TRUE, n = n_nons)
      weights_strap <- weights[strap]
      X_nons_strap <- X_nons[strap,]
      y_strap <- y[strap]

      #using svy package
      strap_rand_svy <- which(rep_weights[,k] != 0)
      weights_rand_strap_svy <- rep_weights[,k] * weights_rand
      N_strap <- sum(weights_rand_strap_svy)

      model_strap <- stats::glm.fit(x = X_nons_strap,
                                    y = y_strap,
                                    weights = weights_strap,
                                    family = family)

      beta <- model_strap$coefficients

      ystrap_rand <- as.numeric(X_rand %*% beta) # TODO with predict.glm

      #mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
      mu_hat_boot <- weighted.mean(x = ystrap_rand, w = weights_rand_strap_svy)
      mu_hats[k] <- mu_hat_boot
      k <- k + 1
    }
  } else if (method == "nn") {

    while (k <= num_boot) {

      strap <- sample.int(replace = TRUE, n = n_nons)
      weights_strap <- weights[strap]
      X_nons_strap <- X_nons[strap,]
      y_strap <- y[strap]

      strap_rand <- sample.int(replace = TRUE, n = n_rand)
      weights_rand_strap <- weights_rand[strap_rand]
      X_rand_strap <- X_rand[strap_rand,]
      N_strap <- sum(weights_rand_strap)

      model_rand <- nonprobMI_nn(data = X_nons_strap,
                                 query = X_rand_strap,
                                 k = control$k,
                                 treetype = control$treetype,
                                 searchtype = control$searchtype)
      y_rand_strap <- vector(mode = "numeric", length = n_rand)

      for (i in 1:n_rand) {
        idx <- model_rand$nn.idx[i,]
        y_rand_strap[i] <- mean(y_strap[idx])
      }

      mu_hat_boot <- weighted.mean(x = ystrap_rand, w = weights_rand_strap_svy)
      mu_hats[k] <- mu_hat_boot
      k <- k + 1
    }

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  boot_var
}

bootIPW <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    R,
                    theta_hat,
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
                    pop_totals = NULL,
                    control_selection,
                    ...){
  mu_hats <- vector(mode = "numeric", length = num_boot)
  if (!is.null(weights_rand)) N <- sum(weights_rand)
  estimation_method <- get_method(est_method)
  method <- get_method(method_selection)
  inv_link <- method$make_link_inv
  k <- 1

  while (k <= num_boot) {

    if (is.null(pop_totals)) {
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
                                      maxit = maxit,
                                      control_selection = control_selection)

      est_method_obj <- estimation_method$estimation_model(model = model_sel,
                                                           method_selection = method_selection)

      ps_nons <- est_method_obj$ps_nons
      weights_nons <- 1/ps_nons
      N_est_nons <- ifelse(is.null(pop_size), sum(weights[strap_nons] * 1/ps_nons), pop_size)

      mu_hat_boot <- mu_hatIPW(y = y[strap_nons],
                               weights = weights[strap_nons],
                               weights_nons = weights_nons,
                               N = N_est_nons) # IPW estimator
      mu_hats[k] <- mu_hat_boot

    } else {
      strap <- sample.int(replace = TRUE, n = n_nons)

      X_strap <- X_nons[strap, ]
      R_strap <- R[strap]
      weights_strap <- weights[strap]

      h_object_strap <- theta_h_estimation(R = R_strap,
                                           X = X_strap,
                                           weights_rand = NULL,
                                           weights = weights_strap,
                                           h = h,
                                           method_selection = method_selection,
                                           maxit = maxit,
                                           pop_totals = pop_totals)
      theta_hat_strap <- h_object_strap$theta_h
      ps_nons <- inv_link(theta_hat_strap %*% t(X_strap))

      weights_nons <- 1/ps_nons
      N_est_nons <- ifelse(is.null(pop_size), sum(weights_strap * weights_nons), pop_size)

      mu_hat_boot <- mu_hatIPW(y = y[strap],
                               weights = weights_strap,
                               weights_nons = weights_nons,
                               N = N_est_nons) # IPW estimator
      mu_hats[k] <- mu_hat_boot
    }
    k <- k + 1
  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  list(boot_var = boot_var)
}

bootDR <- function(SelectionModel,
                   OutcomeModel,
                   family_outcome,
                   num_boot,
                   weights,
                   weights_rand,
                   R,
                   theta_hat,
                   mu_hat,
                   method_selection,
                   control_selection,
                   n_nons,
                   n_rand,
                   optim_method,
                   est_method,
                   h,
                   maxit,
                   pop_size,
                   pop_totals,
                   pop_means,
                   ...) {

  mu_hats <- vector(mode = "numeric", length = num_boot)
  estimation_method <- get_method(est_method)
  k <- 1
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (is.null(pop_totals)) {
    N <- sum(weights_rand)
    while (k <= num_boot) {
        strap_nons <- sample.int(replace = TRUE, n = n_nons)
        strap_rand <- sample.int(replace = TRUE, n = n_rand)

        model_out <- stats::glm.fit(x = OutcomeModel$X_nons[strap_nons, ],
                                    y = OutcomeModel$y[strap_nons],
                                    weights = weights[strap_nons],
                                    family = family)


        model_nons_coefs <- model_out$coefficients
        y_rand_pred <- as.numeric(OutcomeModel$X_rand[strap_rand, ] %*% model_nons_coefs) # TODO with predict.glm
        y_nons_pred <- model_out$fitted.values

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
                                        maxit = maxit,
                                        control_selection = control_selection)

        est_method_obj <- estimation_method$estimation_model(model = model_sel,
                                                             method_selection = method_selection)
        ps_nons <- est_method_obj$ps_nons
        weights_nons <- 1/ps_nons
        N_est_nons <- sum(weights_nons)
        N_est_rand <- sum(weights_rand[strap_rand])

        mu_hat_boot <- mu_hatDR(y = OutcomeModel$y_nons[strap_nons],
                                y_nons = y_nons_pred,
                                y_rand = y_rand_pred,
                                weights = weights[strap_nons],
                                weights_nons = weights_nons,
                                weights_rand = weights_rand[strap_rand],
                                N_nons = N_est_nons,
                                N_rand = N_est_rand)
        mu_hats[k] <- mu_hat_boot
        k <- k + 1
      }
    } else { # TODO
      while (k <= num_boot) {
        #stop("Bootstrap with pop_totals is not yet implemented.")

        strap <- sample.int(replace = TRUE, n = n_nons)
        X_strap <- SelectionModel$X_nons[strap, ]
        y_strap <- SelectionModel$y_nons[strap]
        R_strap <- rep(1, nrow(X_strap))
        weights_strap <- weights[strap]


        h_object_strap <- theta_h_estimation(R = R_strap,
                                             X = X_strap,
                                             weights = weights_strap,
                                             h = h,
                                             method_selection = method_selection,
                                             maxit = maxit,
                                             pop_totals = pop_totals,
                                             weights_rand = NULL)

        theta_hat_strap <- h_object_strap$theta_h
        method <- get_method(method_selection)
        inv_link <- method$make_link_inv
        ps_nons_strap <- inv_link(theta_hat_strap %*% t(X_strap))
        weights_nons_strap <- 1/ps_nons_strap
        N_est <- sum(weights_strap * weights_nons_strap)
        if(is.null(pop_size)) pop_size <- N_est

        model_out_strap <- stats::glm.fit(x = OutcomeModel$X_nons[strap, ],
                                    y = OutcomeModel$y[strap],
                                    weights = weights_strap,
                                    family = family)

        model_nons_coefs <- model_out_strap$coefficients
        y_rand_pred <- as.numeric(pop_totals %*% model_nons_coefs) # TODO with predict.glm
        y_nons_pred <- model_out_strap$fitted.values

        mu_hat_boot <- 1/N_est * sum(weights_nons_strap * (weights_strap * (OutcomeModel$y[strap] - y_nons_pred))) + 1/pop_size * y_rand_pred
        mu_hats[k] <- mu_hat_boot
        k <- k + 1
      }
  }
  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  list(boot_var = boot_var)
}

bootDR_sel <- function(X,
                       R,
                       y,
                       prior_weights,
                       method_selection,
                       family_nonprobsvy,
                       mu_hat,
                       n_nons,
                       n_rand,
                       num_boot,
                       par0,
                       psel) { # TODO function to test
  mu_hats <- vector(mode = "numeric", length = num_boot)
  k <- 1
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- X[loc_nons,]
  X_rand <- X[loc_rand,]
  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]
  method <- get_method(method_selection)
  inv_link <- method$make_link_inv
  while (k <= num_boot) {
    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)
    prior_weights_strap <- c(prior_weights[loc_nons][strap_nons], prior_weights[loc_rand][strap_rand])

    X_strap <- rbind(X_nons[strap_nons, ], X_rand[strap_rand, ])
    y_strap <- c(y_nons[strap_nons], y_rand[strap_rand])


    multiroot <- nleqslv::nleqslv(x = par0, # TODO add user-specified parameters to control functions
                                  fn = u_theta_beta_dr,
                                  method = "Newton", # TODO consider the method
                                  global = "qline",
                                  xscalm = "fixed",
                                  jacobian = TRUE,
                                  control = list(scalex = rep(1, length(par0))), # TODO algorithm did not converge in maxit iterations for cloglog
                                  R = R,
                                  X = X_strap,
                                  y = y_strap,
                                  weights = prior_weights_strap,
                                  method_selection = method_selection,
                                  family_nonprobsvy = family_nonprobsvy)

    X_design_strap <- cbind(1, X_strap)
    par_sel_strap <- multiroot$x

    theta_hat_strap <- par_sel_strap[1:(psel+1)]
    beta_hat_strap <- par_sel_strap[(psel+2):(2*psel+2)]
    ps_strap <- inv_link(as.vector(X_design_strap %*% as.matrix(theta_hat_strap)))
    weights_nons_strap <- 1/ps_strap[loc_nons]
    N_nons <- sum(prior_weights[loc_nons][strap_nons] * weights_nons_strap)
    N_rand <- sum(prior_weights[loc_rand][strap_rand])

    eta <- as.vector(X_design_strap %*% as.matrix(beta_hat_strap))
    y_hat_strap <- family_nonprobsvy$mu(eta)
    y_rand_pred_strap <- y_hat_strap[loc_rand]
    y_nons_pred_strap <- y_hat_strap[loc_nons]

    mu_hat_boot <- mu_hatDR(y = y_nons[strap_nons],
                            y_nons = y_nons_pred_strap,
                            y_rand = y_rand_pred_strap,
                            weights = prior_weights[loc_nons][strap_nons],
                            weights_nons = weights_nons_strap,
                            weights_rand = prior_weights[loc_rand][strap_rand],
                            N_nons = N_nons,
                            N_rand = N_rand) #DR estimator
    mu_hats[k] <- mu_hat_boot
    k <- k + 1
  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  list(boot_var = boot_var)
}
