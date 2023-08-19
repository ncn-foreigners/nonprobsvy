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
                   pop_totals,
                   ...
                   ){ # TODO add methods instead of conditional loops

  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  k <- 1
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
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
        eta <- X_rand %*% beta
        y_strap_rand <- family_nonprobsvy$mu(eta)

        #mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
        mu_hat_boot <- weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
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

        y_rand_strap <- apply(model_rand$nn.idx, 1,
                             FUN=\(x) mean(y_strap[x])
                             #FUN=\(x) mean(sample_nonprob$short_[x])
        )

        mu_hat_boot <- weighted.mean(x = y_rand_strap, w = weights_rand_strap)
        mu_hats[k] <- mu_hat_boot
        k <- k + 1
      }
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      while (k <= num_boot) {
        strap <- sample.int(replace = TRUE, n = n_nons)
        weights_strap <- weights[strap]
        X_nons_strap <- X_nons[strap,]
        y_strap <- y[strap]

        model_strap <- stats::glm.fit(x = X_nons_strap,
                                      y = y_strap,
                                      weights = weights_strap,
                                      family = family)

        beta <- model_strap$coefficients
        eta <- pop_totals %*% beta
        y_strap_rand <- family_nonprobsvy$mu(eta)

        #mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
        mu_hat_boot <-  as.vector(y_strap_rand/N)
        mu_hats[k] <- mu_hat_boot
        k <- k + 1
      }
    } else if (method == "nn") {
      while (k <= num_boot) {
        strap <- sample.int(replace = TRUE, n = n_nons)
        weights_strap <- weights[strap]
        X_nons_strap <- X_nons[strap,]
        y_strap <- y[strap]

        model_rand <- nonprobMI_nn(data = X_nons_strap,
                                   query = t(pop_totals / N),
                                   k = control$k,
                                   treetype = control$treetype,
                                   searchtype = control$searchtype)
        mu_hat_boot <- mean(y_strap[model_rand$nn.idx])
        mu_hats[k] <- mu_hat_boot
        k <- k + 1
      }
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
  k <- 1
  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (est_method == "mm") {
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
    p <- ncol(X)
    y_rand <- vector(mode = "numeric", length = n_rand)
    y <- c(y_rand, OutcomeModel$y_nons) # outcome variable for joint model
    var_obj <- bootDR_sel(X = X,
                          R = R,
                          y = y,
                          prior_weights = c(weights_rand, weights),
                          method_selection = method_selection,
                          family_nonprobsvy = family_nonprobsvy,
                          mu_hat = mu_hat,
                          n_nons = n_nons,
                          n_rand = n_rand,
                          num_boot = num_boot,
                          par0 = rep(0, 2*p),
                          psel = p)
    boot_var <- var_obj$boot_var
  } else {
    estimation_method <- get_method(est_method)
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
          eta <- OutcomeModel$X_rand[strap_rand, ] %*% model_nons_coefs
          y_rand_pred <- family_nonprobsvy$mu(eta)
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
          y#_rand_pred <- as.numeric(pop_totals %*% model_nons_coefs) # TODO with predict.glm
          eta <- OutcomeModel$X_rand[strap_rand, ] %*% model_nons_coefs
          y_rand_pred <- family_nonprobsvy$mu(eta)
          y_nons_pred <- model_out_strap$fitted.values

          mu_hat_boot <- 1/N_est * sum(weights_nons_strap * (weights_strap * (OutcomeModel$y[strap] - y_nons_pred))) + 1/pop_size * y_rand_pred
          mu_hats[k] <- mu_hat_boot
          k <- k + 1
        }
      }
    boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
    }
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
  while (k <= num_boot) {
    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)

    X_strap <- rbind(X_rand[strap_rand, ], X_nons[strap_nons, ])
    y_strap <- c(y_rand[strap_rand], y_nons[strap_nons])

    model_strap <- mm(X = X_strap,
                      y = y_strap,
                      weights = prior_weights[loc_nons][strap_nons],
                      weights_rand = prior_weights[loc_rand][strap_rand],
                      R = R, #c(R[loc_nons][strap_nons], R[loc_rand][strap_rand]),
                      n_nons = n_nons,
                      n_rand = n_rand,
                      method_selection = method_selection,
                      family = family_nonprobsvy,
                      boot = TRUE)

    weights_nons_strap <- 1/model_strap$selection$ps_nons
    N_nons <- sum(prior_weights[loc_nons][strap_nons] * weights_nons_strap)
    N_rand <- sum(prior_weights[loc_rand][strap_rand])

    mu_hat_boot <- mu_hatDR(y = y_nons[strap_nons],
                            y_nons = model_strap$outcome$y_nons_pred,
                            y_rand = model_strap$outcome$y_rand_pred,
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

# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
bootMI_multicore <- function(X_rand,
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
                             pop_totals,
                             cores,
                             ...) {
  #mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
    if (method == "glm") {
      rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "mu_hatIPW", "probit", "cloglog"))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
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
          eta <- X_rand %*% beta
          y_strap_rand <- family_nonprobsvy$mu(eta)

          weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
        })
    } else if (method == "nn") {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "nonprobMI_nn", "probit", "cloglog"))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {

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

          y_rand_strap <- apply(model_rand$nn.idx, 1,
                                FUN=\(x) mean(y_strap[x])
                                #FUN=\(x) mean(sample_nonprob$short_[x])
          )

          weighted.mean(x = y_rand_strap, w = weights_rand_strap)
        })
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "nonprobMI_nn", "probit", "cloglog"))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap,]
          y_strap <- y[strap]

          model_strap <- stats::glm.fit(x = X_nons_strap,
                                        y = y_strap,
                                        weights = weights_strap,
                                        family = family)

          beta <- model_strap$coefficients
          eta <- pop_totals %*% beta
          y_strap_rand <- family_nonprobsvy$mu(eta)

          #mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
          as.vector(y_strap_rand/N)
        })
    } else if (method == "nn") {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "nonprobMI_nn", "probit", "cloglog"))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap,]
          y_strap <- y[strap]

          model_rand <- nonprobMI_nn(data = X_nons_strap,
                                     query = t(pop_totals / N),
                                     k = control$k,
                                     treetype = control$treetype,
                                     searchtype = control$searchtype)
          mean(y_strap[model_rand$nn.idx])
        })
    }
  }
  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  boot_var
}

#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
bootIPW_multicore <- function(X_rand,
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
                              cores,
                              ...) {

  if (!is.null(weights_rand)) N <- sum(weights_rand)
  estimation_method <- get_method(est_method)
  method <- get_method(method_selection)

  inv_link <- method$make_link_inv

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "mu_hatIPW", "probit", "cloglog"))

  mu_hats <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:num_boot, .combine = c),
    ex = {

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


        model_sel <- estimation_method$model_selection(X,
                                                       X_nons[strap_nons, ],
                                                       X_rand[strap_rand, ],
                                                       weights = weights[strap_nons],
                                                       weights_rand = weights_rand[strap_rand],
                                                       R,
                                                       method_selection,
                                                       optim_method,
                                                       h,
                                                       est_method,
                                                       maxit,
                                                       FALSE,
                                                       control_selection)

        est_method_obj <- estimation_method$estimation_model(model = model_sel,
                                                             method_selection = method_selection)

        ps_nons <- est_method_obj$ps_nons
        weights_nons <- 1/ps_nons
        N_est_nons <- ifelse(is.null(pop_size), sum(weights[strap_nons] * 1/ps_nons), pop_size)

        mu_hat_boot <- mu_hatIPW(y = y[strap_nons],
                                 weights = weights[strap_nons],
                                 weights_nons = weights_nons,
                                 N = N_est_nons) # IPW estimator
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
      }

      mu_hat_boot
  })

  # można zmienić na num_boot-1 :)
  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  list(boot_var = boot_var)
}

#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
bootDR_multicore <- function(SelectionModel,
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
                             cores,
                             ...) {

  # mu_hats <- vector(mode = "numeric", length = num_boot)
  # k <- 1
  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (est_method == "mm") {
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
    p <- ncol(X)
    y_rand <- vector(mode = "numeric", length = n_rand)
    y <- c(y_rand, OutcomeModel$y_nons) # outcome variable for joint model
    var_obj <- bootDR_sel_multicore(X = X,
                                    R = R,
                                    y = y,
                                    prior_weights = c(weights_rand, weights),
                                    method_selection = method_selection,
                                    family_nonprobsvy = family_nonprobsvy,
                                    mu_hat = mu_hat,
                                    n_nons = n_nons,
                                    n_rand = n_rand,
                                    num_boot = num_boot,
                                    par0 = rep(0, 2*p),
                                    psel = p,
                                    cores = cores)
    boot_var <- var_obj$boot_var
  } else {
    if (is.null(pop_totals)) {
      N <- sum(weights_rand)

      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "mu_hatDR", "probit", "cloglog"))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          estimation_method <- get_method(est_method)
          strap_nons <- sample.int(replace = TRUE, n = n_nons)
          strap_rand <- sample.int(replace = TRUE, n = n_rand)

          model_out <- stats::glm.fit(x = OutcomeModel$X_nons[strap_nons, ],
                                      y = OutcomeModel$y[strap_nons],
                                      weights = weights[strap_nons],
                                      family = family)


          model_nons_coefs <- model_out$coefficients
          eta <- OutcomeModel$X_rand[strap_rand, ] %*% model_nons_coefs
          y_rand_pred <- family_nonprobsvy$mu(eta)
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

          mu_hatDR(y = OutcomeModel$y_nons[strap_nons],
                   y_nons = y_nons_pred,
                   y_rand = y_rand_pred,
                   weights = weights[strap_nons],
                   weights_nons = weights_nons,
                   weights_rand = weights_rand[strap_rand],
                   N_nons = N_est_nons,
                   N_rand = N_est_rand)
      })
    } else { # TODO
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))

      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
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
        y#_rand_pred <- as.numeric(pop_totals %*% model_nons_coefs) # TODO with predict.glm
        eta <- OutcomeModel$X_rand[strap_rand, ] %*% model_nons_coefs
        y_rand_pred <- family_nonprobsvy$mu(eta)
        y_nons_pred <- model_out_strap$fitted.values

        1/N_est * sum(weights_nons_strap * (weights_strap * (OutcomeModel$y[strap] - y_nons_pred))) + 1/pop_size * y_rand_pred
      })
    }
    boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  }
  list(boot_var = boot_var)
}

# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
bootDR_sel_multicore <- function(X,
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
                                 psel,
                                 cores) { # TODO function to test
  mu_hats <- vector(mode = "numeric", length = num_boot)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- X[loc_nons,]
  X_rand <- X[loc_rand,]
  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl = cl, varlist = c("internal_selection", "logit", "start_fit", "get_method", "controlSel", "mle", "mu_hatIPW", "probit", "cloglog"))

  mu_hats <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:num_boot, .combine = c),
    ex = {
    strap_nons <- sample.int(replace = TRUE, n = n_nons)
    strap_rand <- sample.int(replace = TRUE, n = n_rand)

    X_strap <- rbind(X_rand[strap_rand, ], X_nons[strap_nons, ])
    y_strap <- c(y_rand[strap_rand], y_nons[strap_nons])

    model_strap <- mm(X = X_strap,
                      y = y_strap,
                      weights = prior_weights[loc_nons][strap_nons],
                      weights_rand = prior_weights[loc_rand][strap_rand],
                      R = R, #c(R[loc_nons][strap_nons], R[loc_rand][strap_rand]),
                      n_nons = n_nons,
                      n_rand = n_rand,
                      method_selection = method_selection,
                      family = family_nonprobsvy,
                      boot = TRUE)

    weights_nons_strap <- 1/model_strap$selection$ps_nons
    N_nons <- sum(prior_weights[loc_nons][strap_nons] * weights_nons_strap)
    N_rand <- sum(prior_weights[loc_rand][strap_rand])

    mu_hatDR(y = y_nons[strap_nons],
             y_nons = model_strap$outcome$y_nons_pred,
             y_rand = model_strap$outcome$y_rand_pred,
             weights = prior_weights[loc_nons][strap_nons],
             weights_nons = weights_nons_strap,
             weights_rand = prior_weights[loc_rand][strap_rand],
             N_nons = N_nons,
             N_rand = N_rand) #DR estimator
  })
  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)
  list(boot_var = boot_var)
}
