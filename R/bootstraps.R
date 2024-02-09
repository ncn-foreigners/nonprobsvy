# These functions are only used internally, so there is no need for documenting them
#' @importFrom survey as.svrepdesign
#' @importFrom nleqslv nleqslv
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar

bootMI <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   start_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   svydesign,
                   model_obj = model_obj,
                   rep_type,
                   method,
                   control_outcome,
                   control_inference,
                   pop_totals,
                   verbose,
                   ...) { # TODO add methods instead of conditions

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

  if (is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
  }

<<<<<<< HEAD
  pmm_exact_se <- control_inference$pmm_exact_se
  comp3_stat <- numeric(length = num_boot)
=======
  predictive_match = control_outcome$predictive_match
  pmm_exact_se = control_inference$pmm_exact_se
  pmm_reg_engine = control_outcome$pmm_reg_engine
  pi_ij = control_inference$pi_ij
  pmm_exact_se <- control_inference$pmm_exact_se
  #comp3_stat <- numeric(length = num_boot)
>>>>>>> pmm

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
    rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
    if (method == "glm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            # X_rand_strap <- X_rand[which(rep_weights[,k] != 0),]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta <- X_rand %*% beta
            y_strap_rand <- family_nonprobsvy$linkinv(eta)

            # mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
            mu_hat_boot <- weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "nn") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
            # weights_rand_strap <- weights_rand[strap_rand]
            # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]
            # N_strap <- sum(weights_rand_strap)

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            X_rand_strap <- X_rand[which(rep_weights[, k] != 0), ]
            weights_rand_strap <- weights_rand_strap_svy[strap_rand_svy]

            model_rand <- nonprobMI_nn(
              data = X_nons_strap,
              query = X_rand_strap,
              k = control_outcome$k,
              treetype = control_outcome$treetype,
              searchtype = control_outcome$searchtype
            )

            y_rand_strap <- apply(model_rand$nn.idx, 1,
              FUN = \(x) mean(y_strap[x])
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )

            mu_hat_boot <- weighted.mean(x = y_rand_strap, w = weights_rand_strap)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "pmm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
            # weights_rand_strap <- weights_rand[strap_rand]
            # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]
            # N_strap <- sum(weights_rand_strap)

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            X_rand_strap <- X_rand[which(rep_weights[, k] != 0), ]
            n_rand_strap <- nrow(X_rand_strap)
            weights_rand_strap <- weights_rand_strap_svy[strap_rand_svy]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta_rand <- X_rand_strap %*% beta
            eta_nons <- X_nons_strap %*% beta
            y_rand_strap <- family_nonprobsvy$linkinv(eta_rand)
            y_nons_strap <- family_nonprobsvy$linkinv(eta_nons)


            model_rand <- switch(control_outcome$predictive_match,
              { # 1
                nonprobMI_nn(
                  data = y_strap,
                  query = y_rand_strap,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              },
              { # 2
                nonprobMI_nn(
                  data = y_nons_strap,
                  query = y_rand_strap,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              }
            )

            y_rand_strap <- apply(model_rand$nn.idx, 1,
              FUN = \(x) mean(y_strap[x])
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )

            mu_hat_boot <- weighted.mean(x = y_rand_strap, w = weights_rand_strap)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }

<<<<<<< HEAD
            comp3 <- 0
            if (pmm_exact_se) {
              pi_ij <- outer(1 / weights_rand_strap, 1 / weights_rand_strap) * (
                1 - outer(1 - 1 / weights_rand_strap, 1 - 1 / weights_rand_strap) / sum(1 - 1 / weights_rand_strap)
              )
              mat_preds <- matrix(y_strap[model_rand$nn.idx[1:n_rand_strap, ]], nrow = n_rand_strap) / control_outcome$k
              for (ii in 1:n_rand_strap) {
                for (jj in 1:ii) {
                  comp3 <- comp3 +
                    ((pi_ij[ii, jj] ^ -1) / N ^ 2) *
                    (sum(outer(mat_preds[ii,], mat_preds[jj, ])) -
                       y_rand_strap[ii] * y_rand_strap[jj])
                }
              }
            }
            comp3_stat[k] <- comp3
=======
            # slower option
            # if (pmm_exact_se) {
            #   comp2 <- pmm_exact(pi_ij,
            #                      weights_rand_strap,
            #                      n_nons = n_nons,
            #                      y = y,
            #                      pmm_reg_engine = pmm_reg_engine,
            #                      model_obj = model_obj,
            #                      svydesign = svydesign,
            #                      predictive_match = predictive_match,
            #                      k = control_inference$k,
            #                      N = N)
            # }

            #comp2_stat[k] <- comp2
>>>>>>> pmm
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta <- pop_totals %*% beta / N
            y_strap_rand <- family_nonprobsvy$linkinv(eta)

            # mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
            mu_hat_boot <- as.vector(y_strap_rand)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "nn") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_rand <- nonprobMI_nn(
              data = X_nons_strap,
              query = t(pop_totals / N),
              k = control_outcome$k,
              treetype = control_outcome$treetype,
              searchtype = control_outcome$searchtype
            )
            mu_hat_boot <- mean(y_strap[model_rand$nn.idx])
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "pmm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta_rand <- pop_totals %*% beta / N
            eta_nons <- X_nons_strap %*% beta
            y_strap_rand <- family_nonprobsvy$linkinv(eta_rand)
            y_strap_nons <- family_nonprobsvy$linkinv(eta_nons)


            model_rand <- switch(control_outcome$predictive_match,
                                 { # 1
                                   nonprobMI_nn(
                                     data = y_strap,
                                     query = y_strap_rand,
                                     k = control_outcome$k,
                                     treetype = control_outcome$treetype,
                                     searchtype = control_outcome$searchtype
                                   )
                                 },
                                 { # 2
                                   nonprobMI_nn(
                                     data = y_strap_nons,
                                     query = y_strap_rand,
                                     k = control_outcome$k,
                                     treetype = control_outcome$treetype,
                                     searchtype = control_outcome$searchtype
                                   )
                                 }
            )
            #
            # model_rand <- nonprobMI_nn(
            #   data = y_strap_nons,
            #   query = y_strap_rand,
            #   k = control_outcome$k,
            #   treetype = control_outcome$treetype,
            #   searchtype = control_outcome$searchtype
            # )


            mu_hat_boot <- mean(y_strap[model_rand$nn.idx])
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    }
  }
  # mu_hat_boot <- mean(mu_hats)
  if (method == "pmm") {
<<<<<<< HEAD
    comp3_mean <- mean(comp3)
  } else {
    comp3_mean <- 0
  }
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2) + comp3_mean
=======
    if (pmm_exact_se) {
      comp2 <- pmm_exact(pi_ij,
                         weights_rand,
                         n_nons = n_nons,
                         y = y,
                         pmm_reg_engine = pmm_reg_engine,
                         model_obj = model_obj,
                         svydesign = svydesign,
                         predictive_match = predictive_match,
                         k = control_inference$k,
                         N = N)
    } else {
      comp2 <- 0
    }
  } else {
    comp2 <- 0
  }
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2) + comp2
>>>>>>> pmm
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats,
<<<<<<< HEAD
    comp3_stat = comp3_stat
=======
    comp2 = comp2
>>>>>>> pmm
  )
}

bootIPW <- function(X_rand,
                    X_nons,
                    svydesign,
                    weights,
                    ys,
                    R,
                    theta_hat,
                    num_boot,
                    weights_rand,
                    mu_hats,
                    method_selection,
                    start_selection,
                    n_nons,
                    n_rand,
                    optim_method,
                    est_method,
                    h,
                    rep_type,
                    maxit,
                    control_inference,
                    control_selection,
                    verbose,
                    pop_size,
                    pop_totals,
                    ...) {
  if (!is.null(weights_rand)) N <- sum(weights_rand)
  estimation_method <- get_method(est_method)
  method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection_function)
  inv_link <- method$make_link_inv
  k <- 1
  rep_type <- control_inference$rep_type
  mu_len <- length(mu_hats)
  mu_hats_boot <- matrix(nrow = num_boot, ncol = mu_len)
  boot_vars <- numeric(length = mu_len)

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
  }

  if (is.null(pop_totals)) {
    rep_weights <- survey::as.svrepdesign(design = svydesign, type = rep_type, replicates = num_boot)$repweights$weights # TODO customise to calibrated svydesign
    while (k <= num_boot) {
      tryCatch(
        {
          strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

          # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
          # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]

          X_nons_strap <- X_nons[strap_nons, , drop = FALSE]
          X <- rbind(X_rand_strap, X_nons_strap)
          n_rand_strap <- nrow(X_rand_strap)

          R_nons <- rep(1, n_nons)
          R_rand <- rep(0, n_rand_strap)
          R <- c(R_rand, R_nons)

          model_sel <- internal_selection(
            X = X,
            X_nons = X_nons_strap,
            X_rand = X_rand_strap,
            weights = weights[strap_nons],
            weights_rand = weights_strap_rand,
            R = R,
            method_selection = method_selection,
            optim_method = optim_method,
            h = h,
            est_method = est_method,
            maxit = maxit,
            control_selection = control_selection,
            start = start_selection
          )

          est_method_obj <- estimation_method$estimation_model(
            model = model_sel,
            method_selection = method_selection
          )

          ps_nons <- est_method_obj$ps_nons
          weights_nons <- 1 / ps_nons
          N_est_nons <- ifelse(is.null(pop_size), sum(weights[strap_nons] * weights_nons), pop_size)

          for (l in 1:mu_len) {
            mu_hats_boot[k, l] <- mu_hatIPW(
              y = ys[[l]][strap_nons],
              weights = weights[strap_nons],
              weights_nons = weights_nons,
              N = N_est_nons
            ) # IPW estimator
          }
          if (verbose) {
            # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hats_boot[k,], sep = "")
            # print(info)
            utils::setTxtProgressBar(pb, k)
          }
          k <- k + 1
        },
        error = function(e) {
          if (verbose) {
            info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
            print(info)
          }
        }
      )
    }
  } else {
    while (k <= num_boot) {
      tryCatch(
        {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)

          X_strap <- X_nons[strap, , drop = FALSE]
          R_strap <- R[strap]
          weights_strap <- weights[strap]

          h_object_strap <- theta_h_estimation(
            R = R_strap,
            X = X_strap,
            weights_rand = NULL,
            weights = weights_strap,
            h = h,
            method_selection = method_selection,
            maxit = maxit,
            pop_totals = pop_totals,
            start = start_selection
          )
          theta_hat_strap <- h_object_strap$theta_h
          ps_nons <- inv_link(theta_hat_strap %*% t(X_strap))

          weights_nons <- 1 / ps_nons
          N_est_nons <- ifelse(is.null(pop_size), sum(weights_strap * weights_nons), pop_size)

          for (l in 1:mu_len) {
            mu_hats_boot[k, l] <- mu_hatIPW(
              y = ys[[l]][strap],
              weights = weights_strap,
              weights_nons = weights_nons,
              N = N_est_nons
            ) # IPW estimator
          }
          if (verbose) {
            # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hats_boot[k], sep = "")
            # print(info)
            utils::setTxtProgressBar(pb, k)
          }
          k <- k + 1
        },
        error = function(e) {
          if (verbose) {
            info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
            print(info)
          }
        }
      )
    }
  }
  # mu_hats_boot_means <- colMeans(mu_hats_boot)
  # boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat_boot)^2)
  for (l in 1:mu_len) {
    boot_vars[l] <- 1 / (num_boot - 1) * sum((mu_hats_boot[, l] - mu_hats[l])^2)
  }
  if (verbose) {
    close(pb)
  }
  list(
    var = boot_vars,
    # mu = mu_hats_boot_means,
    stat = mu_hats_boot
  )
}

bootDR <- function(outcome,
                   data,
                   svydesign,
                   SelectionModel,
                   OutcomeModel,
                   family_outcome,
                   method_outcome,
                   start_outcome,
                   num_boot,
                   weights,
                   weights_rand,
                   R,
                   theta_hat,
                   mu_hat,
                   method_selection,
                   start_selection,
                   control_selection,
                   control_outcome,
                   control_inference,
                   n_nons,
                   n_rand,
                   optim_method,
                   est_method,
                   h,
                   maxit,
                   pop_size,
                   pop_totals,
                   pop_means,
                   bias_correction,
                   verbose,
                   ...) {
  mu_hats <- vector(mode = "numeric", length = num_boot)
  k <- 1
  rep_type <- control_inference$rep_type
  if (is.character(family_outcome)) {
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
  method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
  MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())

  if (bias_correction == TRUE) {
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
    p <- ncol(X)
    y_rand <- vector(mode = "numeric", length = n_rand)
    y <- c(y_rand, OutcomeModel$y_nons) # outcome variable for joint model
    var_obj <- bootDR_sel(
      X = X,
      R = R,
      y = y,
      svydesign = svydesign,
      rep_type = rep_type,
      weights = weights,
      weights_rand = weights_rand,
      method_selection = method_selection,
      family_nonprobsvy = family_nonprobsvy,
      mu_hat = mu_hat,
      n_nons = n_nons,
      n_rand = n_rand,
      num_boot = num_boot,
      start_selection = start_selection,
      start_outcome = start_outcome,
      verbose = verbose
    )
    boot_var <- var_obj$var
    mu_hat_boot <- var_obj$mu
  } else {
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
    }
    estimation_method <- get_method(est_method)
    if (is.null(pop_totals)) {
      rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
      N <- sum(weights_rand)
      while (k <= num_boot) {
        tryCatch(
          {
            strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            # N_strap <- sum(weights_rand_strap_svy)
            # X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
            weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

            model_obj <- MethodOutcome(
              outcome = outcome,
              data = data[strap_nons, ],
              weights = weights[strap_nons],
              family_outcome = family_outcome,
              start_outcome = start_outcome,
              X_nons = OutcomeModel$X_nons[strap_nons, , drop = FALSE],
              y_nons = OutcomeModel$y_nons[strap_nons],
              X_rand = OutcomeModel$X_rand[strap_rand_svy, , drop = FALSE],
              control = control_outcome,
              n_nons = n_nons,
              n_rand = n_rand,
              model_frame = OutcomeModel$model_frame_rand[strap_rand_svy, ],
              vars_selection = control_inference$vars_selection,
              pop_totals = pop_totals
            )


            y_rand_pred <- model_obj$y_rand_pred
            y_nons_pred <- model_obj$y_nons_pred

            X_sel <- rbind(
              SelectionModel$X_rand[strap_rand_svy, , drop = FALSE],
              SelectionModel$X_nons[strap_nons, , drop = FALSE]
            )
            n_rand_strap <- nrow(SelectionModel$X_rand[strap_rand_svy, , drop = FALSE])

            R_nons <- rep(1, n_nons)
            R_rand <- rep(0, n_rand_strap)
            R <- c(R_rand, R_nons)

            model_sel <- internal_selection(
              X = X_sel,
              X_nons = SelectionModel$X_nons[strap_nons, , drop = FALSE],
              X_rand = SelectionModel$X_rand[strap_rand_svy, , drop = FALSE],
              weights = weights[strap_nons],
              weights_rand = weights_strap_rand,
              R = R,
              method_selection = method_selection,
              optim_method = optim_method,
              h = h,
              est_method = est_method,
              maxit = maxit,
              control_selection = control_selection,
              start = start_selection
            )

            est_method_obj <- estimation_method$estimation_model(
              model = model_sel,
              method_selection = method_selection
            )
            ps_nons <- est_method_obj$ps_nons
            weights_nons <- 1 / ps_nons
            N_est_nons <- sum(weights_nons)
            N_est_rand <- sum(weights_strap_rand)

            mu_hat_boot <- mu_hatDR(
              y = OutcomeModel$y_nons[strap_nons],
              y_nons = y_nons_pred,
              y_rand = y_rand_pred,
              weights = weights[strap_nons],
              weights_nons = weights_nons,
              weights_rand = weights_strap_rand,
              N_nons = N_est_nons,
              N_rand = N_est_rand
            )
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            X_strap_nons <- SelectionModel$X_nons[strap, , drop = FALSE]
            y_strap <- OutcomeModel$y_nons[strap]
            R_strap <- rep(1, n_nons)
            weights_strap <- weights[strap]
            n_rand <- 0
            X_strap_rand <- NULL

            h_object_strap <- theta_h_estimation(
              R = R_strap,
              X = X_strap_nons,
              weights = weights_strap,
              h = h,
              method_selection = method_selection,
              maxit = maxit,
              pop_totals = pop_totals,
              weights_rand = NULL,
              start = start_selection
            )

            theta_hat_strap <- h_object_strap$theta_h
            method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
            method <- get_method(method_selection_function)
            inv_link <- method$make_link_inv
            ps_nons_strap <- inv_link(theta_hat_strap %*% t(X_strap_nons))
            weights_nons_strap <- 1 / ps_nons_strap
            N_est <- sum(weights_strap * weights_nons_strap)
            if (is.null(pop_size)) pop_size <- N_est

            model_obj <- MethodOutcome(
              outcome = outcome,
              data = data[strap, , drop = FALSE],
              weights = weights_strap,
              family_outcome = family_outcome,
              start_outcome = start_outcome,
              X_nons = X_strap_nons,
              y_nons = y_strap,
              X_rand = X_strap_rand,
              control = control_outcome,
              n_nons = n_nons,
              n_rand = n_rand,
              model_frame = OutcomeModel$model_frame_rand,
              vars_selection = control_inference$vars_selection,
              pop_totals = pop_totals
            )

            y_rand_pred <- model_obj$y_rand_pred
            y_nons_pred <- model_obj$y_nons_pred

            mu_hat_boot <- 1 / N_est * sum(weights_nons_strap * (weights_strap * (y_strap - y_nons_pred))) + ifelse(method_outcome == "glm", 1 / pop_size * y_rand_pred, y_rand_pred)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    }
    # mu_hat_boot <- mean(mu_hats)
    boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  }
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}

bootDR_sel <- function(X,
                       R,
                       y,
                       svydesign,
                       weights,
                       weights_rand,
                       method_selection,
                       family_nonprobsvy,
                       mu_hat,
                       n_nons,
                       n_rand,
                       num_boot,
                       rep_type,
                       start_selection,
                       start_outcome,
                       verbose) { # TODO function to test
  mu_hats <- vector(mode = "numeric", length = num_boot)
  k <- 1
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- X[loc_nons, , drop = FALSE]
  X_rand <- X[loc_rand, , drop = FALSE]
  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
  }

  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
  while (k <= num_boot) {
    tryCatch(
      {
        strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
        # strap_rand <- sample.int(replace = TRUE, n = n_rand,  prob = 1/weights_rand)

        weights_nons <- weights[strap_nons]
        # weights_rand_strap <- weights_rand[strap_rand]

        # using svy package
        strap_rand_svy <- which(rep_weights[, k] != 0)
        weights_rand_strap_svy <- rep_weights[, k] * weights_rand
        N_strap <- sum(weights_rand_strap_svy)
        # X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
        weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]


        X_strap <- rbind(X_rand[strap_rand_svy, , drop = FALSE], X_nons[strap_nons, , drop = FALSE])
        y_strap <- c(y_rand[strap_rand_svy], y_nons[strap_nons])
        n_rand_strap <- nrow(X_rand[strap_rand_svy, , drop = FALSE])

        R_nons <- rep(1, n_nons)
        R_rand <- rep(0, n_rand_strap)
        R <- c(R_rand, R_nons)

        model_strap <- mm(
          X = X_strap,
          y = y_strap,
          weights = weights_nons,
          weights_rand = weights_strap_rand,
          R = R, # c(R[loc_nons][strap_nons], R[loc_rand][strap_rand]),
          n_nons = n_nons,
          n_rand = n_rand_strap,
          method_selection = method_selection,
          family = family_nonprobsvy,
          start_outcome = start_outcome,
          start_selection = start_selection,
          boot = TRUE
        )

        weights_nons_strap <- 1 / model_strap$selection$ps_nons
        N_nons <- sum(weights_nons * weights_nons_strap)
        N_rand <- sum(weights_strap_rand)

        mu_hat_boot <- mu_hatDR(
          y = y_nons[strap_nons],
          y_nons = model_strap$outcome$y_nons_pred,
          y_rand = model_strap$outcome$y_rand_pred,
          weights = weights_nons,
          weights_nons = weights_nons_strap,
          weights_rand = weights_strap_rand,
          N_nons = N_nons,
          N_rand = N_rand
        ) # DR estimator
        mu_hats[k] <- mu_hat_boot
        if (verbose) {
          # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
          # print(info)
          utils::setTxtProgressBar(pb, k)
        }
        k <- k + 1
      },
      error = function(e) {
        if (verbose) {
          info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
          print(info)
        }
      }
    )
  }
  # mu_hat_boot <- mean(mu_hats)
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}

# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
bootMI_multicore <- function(X_rand,
                             X_nons,
                             weights,
                             y,
                             family_outcome,
                             start_outcome,
                             num_boot,
                             weights_rand,
                             mu_hat,
                             svydesign,
                             method,
                             control_outcome,
                             control_inference,
                             pop_totals,
                             cores,
                             verbose,
                             ...) {
  # mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  rep_type <- control_inference$rep_type

  if (is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))

  ## progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(total = num_boot)
    opts <- list(progress = \(n) pb$tick())
  } else {
    opts <- NULL
  }
  ###
  parallel::clusterExport(cl = cl, varlist = c(
    "internal_selection", "logit_model_nonprobsvy", "start_fit", "get_method", "controlSel",
    "mle", "mu_hatIPW", "probit_model_nonprobsvy", "cloglog_model_nonprobsvy", "nonprobMI_nn"
  ))

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
    rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
    if (method == "glm") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          # X_rand_strap <- X_rand[which(rep_weights[,k] != 0),]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta <- X_rand %*% beta
          y_strap_rand <- family_nonprobsvy$linkinv(eta)
          weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
        }
      )
    } else if (method == "nn") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
          # weights_rand_strap <- weights_rand[strap_rand]
          # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]
          # N_strap <- sum(weights_rand_strap)

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

          model_rand <- nonprobMI_nn(
            data = X_nons_strap,
            query = X_rand_strap,
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          y_rand_strap <- apply(model_rand$nn.idx, 1,
            FUN = \(x) mean(y_strap[x])
            # FUN=\(x) mean(sample_nonprob$short_[x])
          )
          weighted.mean(x = y_rand_strap, w = weights_strap_rand)
        }
      )
    } else if (method == "pmm") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
          # weights_rand_strap <- weights_rand[strap_rand]
          # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]
          # N_strap <- sum(weights_rand_strap)

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[weights_rand_strap_svy != 0]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )


          beta <- model_strap$coefficients
          eta_rand <- X_rand_strap %*% beta
          eta_nons <- X_nons_strap %*% beta
          y_rand_strap <- family_nonprobsvy$linkinv(eta_rand)
          y_nons_strap <- family_nonprobsvy$linkinv(eta_nons)


          model_rand <- switch(control_outcome$predictive_match,
            { # 1
              nonprobMI_nn(
                data = y_strap,
                query = y_rand_strap,
                k = control_outcome$k,
                treetype = control_outcome$treetype,
                searchtype = control_outcome$searchtype
              )
            },
            { # 2
              nonprobMI_nn(
                data = y_nons_strap,
                query = y_rand_strap,
                k = control_outcome$k,
                treetype = control_outcome$treetype,
                searchtype = control_outcome$searchtype
              )
            }
          )

          y_rand_strap <- apply(model_rand$nn.idx, 1,
            FUN = \(x) mean(y_strap[x])
            # FUN=\(x) mean(sample_nonprob$short_[x])
          )
          weighted.mean(x = y_rand_strap, w = weights_strap_rand)
        }
      )
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta <- pop_totals %*% beta / N
          y_strap_rand <- family_nonprobsvy$linkinv(eta)

          # mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
          as.vector(y_strap_rand)
        }
      )
    } else if (method == "nn") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_rand <- nonprobMI_nn(
            data = X_nons_strap,
            query = t(pop_totals / N),
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          mean(y_strap[model_rand$nn.idx])
        }
      )
    } else if (method == "pmm") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta_rand <- pop_totals %*% beta
          eta_nons <- X_nons_strap %*% beta
          y_strap_rand <- family_nonprobsvy$linkinv(eta_rand)
          y_strap_nons <- family_nonprobsvy$linkinv(eta_nons)


          model_rand <- nonprobMI_nn(
            data = y_strap_nons,
            query = y_strap_rand,
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          mean(y_strap[model_rand$nn.idx])
        }
      )
    }
  }
  # mu_hat_boot <- mean(mu_hats)
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}

#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
bootIPW_multicore <- function(X_rand,
                              X_nons,
                              svydesign,
                              weights,
                              ys,
                              R,
                              theta_hat,
                              num_boot,
                              weights_rand,
                              mu_hats,
                              method_selection,
                              start_selection,
                              n_nons,
                              n_rand,
                              optim_method,
                              est_method,
                              h,
                              maxit,
                              control_selection,
                              control_inference,
                              cores,
                              verbose,
                              pop_size,
                              pop_totals,
                              ...) {
  if (!is.null(weights_rand)) N <- sum(weights_rand)
  estimation_method <- get_method(est_method)
  method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection_function)
  inv_link <- method$make_link_inv
  rep_type <- control_inference$rep_type

  mu_len <- length(mu_hats)
  mu_hats_boot <- numeric(length = num_boot * mu_len)
  boot_vars <- numeric(length = mu_len)

  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  ## progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(total = num_boot)
    opts <- list(progress = \(n) pb$tick())
  } else {
    opts <- NULL
  }
  ###
  parallel::clusterExport(cl = cl, varlist = c(
    "internal_selection", "logit_model_nonprobsvy", "start_fit", "get_method", "controlSel",
    "mle", "mu_hatIPW", "probit_model_nonprobsvy", "cloglog_model_nonprobsvy", "theta_h_estimation"
  ))

  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

  k <- 1:num_boot
  mu_hats_boot <- foreach::`%dopar%`(
    obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
    ex = {
      if (is.null(pop_totals)) {
        strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)

        # using svy package
        strap_rand_svy <- which(rep_weights[, k] != 0)
        weights_rand_strap_svy <- rep_weights[, k] * weights_rand
        N_strap <- sum(weights_rand_strap_svy)
        X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
        weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

        # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
        # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]

        X_nons_strap <- X_nons[strap_nons, , drop = FALSE]
        X <- rbind(X_rand_strap, X_nons_strap)
        n_rand_strap <- nrow(X_rand_strap)

        R_nons <- rep(1, n_nons)
        R_rand <- rep(0, n_rand_strap)
        R <- c(R_rand, R_nons)

        model_sel <- internal_selection(
          X = X,
          X_nons = X_nons_strap,
          X_rand = X_rand_strap,
          weights = weights[strap_nons],
          weights_rand = weights_strap_rand,
          R = R,
          method_selection = method_selection,
          optim_method = optim_method,
          h = h,
          est_method = est_method,
          maxit = maxit,
          control_selection = control_selection,
          start = start_selection
        )

        est_method_obj <- estimation_method$estimation_model(
          model = model_sel,
          method_selection = method_selection
        )

        ps_nons <- est_method_obj$ps_nons
        weights_nons <- 1 / ps_nons
        N_est_nons <- ifelse(is.null(pop_size), sum(weights[strap_nons] * weights_nons), pop_size)

        # mu_hat_boot <- mu_hatIPW(
        #   y = y[strap_nons],
        #   weights = weights[strap_nons],
        #   weights_nons = weights_nons,
        #   N = N_est_nons
        # ) # IPW estimator

        mu_hats_this_boot <- numeric(mu_len)

        for (l in 1:mu_len) {
          mu_hats_this_boot[l] <- mu_hatIPW(
            y = ys[[l]][strap_nons],
            weights = weights[strap_nons],
            weights_nons = weights_nons,
            N = N_est_nons
          ) # IPW estimator
        }
        mu_hats_this_boot
      } else {
        strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
        X_strap <- X_nons[strap, , drop = FALSE]
        R_strap <- R[strap]
        weights_strap <- weights[strap]

        h_object_strap <- theta_h_estimation(
          R = R_strap,
          X = X_strap,
          weights_rand = NULL,
          weights = weights_strap,
          h = h,
          method_selection = method_selection,
          maxit = maxit,
          pop_totals = pop_totals,
          start = start_selection
        )
        theta_hat_strap <- h_object_strap$theta_h
        ps_nons <- inv_link(theta_hat_strap %*% t(X_strap))

        weights_nons <- 1 / ps_nons
        N_est_nons <- ifelse(is.null(pop_size), sum(weights_strap * weights_nons), pop_size)

        # mu_hat_boot <- mu_hatIPW(
        #   y = y[strap],
        #   weights = weights_strap,
        #   weights_nons = weights_nons,
        #   N = N_est_nons
        # ) # IPW estimator
        for (l in 1:mu_len) {
          mu_hats_boot[k, l] <- mu_hatIPW(
            y = ys[[l]][strap],
            weights = weights_strap,
            weights_nons = weights_nons,
            N = N_est_nons
          ) # IPW estimator
        }
        mu_hats_boot
      }
    }
  )
  mu_hats_boot <- matrix(mu_hats_boot, nrow = num_boot, ncol = mu_len, byrow = TRUE)
  # mu_hats_boot_means <- colMeans(mu_hats_boot)
  # boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat_boot)^2)
  for (l in 1:mu_len) {
    boot_vars[l] <- 1 / (num_boot - 1) * sum((mu_hats_boot[, l] - mu_hats[l])^2)
  }
  list(
    var = boot_vars,
    # mu = mu_hats_boot_means,
    stat = mu_hats_boot
  )
}

#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
bootDR_multicore <- function(outcome,
                             data,
                             svydesign,
                             SelectionModel,
                             OutcomeModel,
                             family_outcome,
                             method_outcome,
                             start_outcome,
                             num_boot,
                             weights,
                             weights_rand,
                             R,
                             theta_hat,
                             mu_hat,
                             method_selection,
                             control_selection,
                             start_selection,
                             control_outcome,
                             control_inference,
                             n_nons,
                             n_rand,
                             optim_method,
                             est_method,
                             h,
                             maxit,
                             pop_size,
                             pop_totals,
                             pop_means,
                             bias_correction,
                             cores,
                             verbose,
                             ...) {
  # mu_hats <- vector(mode = "numeric", length = num_boot)
  # k <- 1
  if (is.character(family_outcome)) {
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
  rep_type <- control_inference$rep_type
  method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
  MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())

  if (bias_correction == TRUE) {
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
    p <- ncol(X)
    y_rand <- vector(mode = "numeric", length = n_rand)
    y <- c(y_rand, OutcomeModel$y_nons) # outcome variable for joint model
    var_obj <- bootDR_sel_multicore(
      X = X,
      R = R,
      y = y,
      svydesign = svydesign,
      rep_type = rep_type,
      weights = weights,
      weights_rand = weights_rand,
      method_selection = method_selection,
      family_nonprobsvy = family_nonprobsvy,
      mu_hat = mu_hat,
      n_nons = n_nons,
      n_rand = n_rand,
      num_boot = num_boot,
      start_selection = start_selection,
      start_outcome = start_outcome,
      cores = cores
    )
    boot_var <- var_obj$var
    mu_hat_boot <- var_obj$mu
  } else {
    rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl))
    ## progress bar
    if (verbose) {
      pb <- progress::progress_bar$new(total = num_boot)
      opts <- list(progress = \(n) pb$tick())
    } else {
      opts <- NULL
    }
    parallel::clusterExport(cl = cl, varlist = c(
      "internal_selection", "internal_outcome", "logit_model_nonprobsvy", "start_fit", "get_method", "controlSel", "theta_h_estimation",
      "mle", "mu_hatDR", "probit_model_nonprobsvy", "cloglog_model_nonprobsvy", "glm_nonprobsvy", "nn_nonprobsvy", "pmm_nonprobsvy",
      "gaussian_nonprobsvy", "poisson_nonprobsvy", "binomial_nonprobsvy", "nonprobMI_fit", "controlOut"
    ))
    if (is.null(pop_totals)) {
      N <- sum(weights_rand)

      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
        ex = {
          estimation_method <- get_method(est_method)
          strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          # N_strap <- sum(weights_rand_strap_svy)
          # X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

          model_obj <- MethodOutcome(
            outcome = outcome,
            data = data[strap_nons, ],
            weights = weights[strap_nons],
            family_outcome = family_outcome,
            start_outcome = start_outcome,
            X_nons = OutcomeModel$X_nons[strap_nons, , drop = FALSE],
            y_nons = OutcomeModel$y_nons[strap_nons],
            X_rand = OutcomeModel$X_rand[strap_rand_svy, , drop = FALSE],
            control = control_outcome,
            n_nons = n_nons,
            n_rand = n_rand,
            model_frame = OutcomeModel$model_frame_rand[strap_rand_svy, ],
            vars_selection = control_inference$vars_selection,
            pop_totals = pop_totals
          )


          y_rand_pred <- model_obj$y_rand_pred
          y_nons_pred <- model_obj$y_nons_pred

          X_sel <- rbind(
            SelectionModel$X_rand[strap_rand_svy, , drop = FALSE],
            SelectionModel$X_nons[strap_nons, , drop = FALSE]
          )
          n_rand_strap <- nrow(SelectionModel$X_rand[strap_rand_svy, , drop = FALSE])

          R_nons <- rep(1, n_nons)
          R_rand <- rep(0, n_rand_strap)
          R <- c(R_rand, R_nons)

          model_sel <- internal_selection(
            X = X_sel,
            X_nons = SelectionModel$X_nons[strap_nons, , drop = FALSE],
            X_rand = SelectionModel$X_rand[strap_rand_svy, , drop = FALSE],
            weights = weights[strap_nons],
            weights_rand = weights_strap_rand,
            R = R,
            method_selection = method_selection,
            optim_method = optim_method,
            h = h,
            est_method = est_method,
            maxit = maxit,
            control_selection = control_selection,
            start = start_selection
          )

          est_method_obj <- estimation_method$estimation_model(
            model = model_sel,
            method_selection = method_selection
          )
          ps_nons <- est_method_obj$ps_nons
          weights_nons <- 1 / ps_nons
          N_est_nons <- sum(weights_nons)
          N_est_rand <- sum(weights_strap_rand)

          mu_hatDR(
            y = OutcomeModel$y_nons[strap_nons],
            y_nons = y_nons_pred,
            y_rand = y_rand_pred,
            weights = weights[strap_nons],
            weights_nons = weights_nons,
            weights_rand = weights_strap_rand,
            N_nons = N_est_nons,
            N_rand = N_est_rand
          )
        }
      )
    } else {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          X_nons_strap <- SelectionModel$X_nons[strap, , drop = FALSE]
          y_strap <- OutcomeModel$y_nons[strap]
          R_strap <- rep(1, n_nons)
          weights_strap <- weights[strap]
          X_rand_strap <- NULL

          h_object_strap <- theta_h_estimation(
            R = R_strap,
            X = X_nons_strap,
            weights = weights_strap,
            h = h,
            method_selection = method_selection,
            maxit = maxit,
            pop_totals = pop_totals,
            start = start_selection,
            weights_rand = NULL
          )

          theta_hat_strap <- h_object_strap$theta_h
          method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
          method <- get_method(method_selection_function)
          inv_link <- method$make_link_inv
          ps_nons_strap <- inv_link(theta_hat_strap %*% t(X_nons_strap))
          weights_nons_strap <- 1 / ps_nons_strap
          N_est <- sum(weights_strap * weights_nons_strap)
          if (is.null(pop_size)) pop_size <- N_est

          model_obj <- MethodOutcome(
            outcome = outcome,
            data = data[strap, , drop = FALSE],
            weights = weights_strap,
            family_outcome = family_outcome,
            start_outcome = start_outcome,
            X_nons = X_nons_strap,
            y_nons = y_strap,
            X_rand = X_rand_strap,
            control = control_outcome,
            n_nons = n_nons,
            n_rand = n_rand,
            model_frame = OutcomeModel$model_frame_rand,
            vars_selection = control_inference$vars_selection,
            pop_totals = pop_totals
          )

          y_rand_pred <- model_obj$y_rand_pred
          y_nons_pred <- model_obj$y_nons_pred

          mu_hat_boot <- 1 / N_est * sum(weights_nons_strap * (weights_strap * (y_strap - y_nons_pred))) + ifelse(method_outcome == "glm", 1 / pop_size * y_rand_pred, y_rand_pred)
          mu_hat_boot
        }
      )
    }
    # mu_hat_boot <- mean(mu_hats)
    boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  }
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}

# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
bootDR_sel_multicore <- function(X,
                                 svydesign,
                                 R,
                                 y,
                                 weights,
                                 weights_rand,
                                 method_selection,
                                 family_nonprobsvy,
                                 mu_hat,
                                 n_nons,
                                 n_rand,
                                 num_boot,
                                 rep_type,
                                 start_selection,
                                 start_outcome,
                                 cores,
                                 verbose) { # TODO function to test
  mu_hats <- vector(mode = "numeric", length = num_boot)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- X[loc_nons, , drop = FALSE]
  X_rand <- X[loc_rand, , drop = FALSE]
  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]

  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  ## progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(total = num_boot)
    opts <- list(progress = \(n) pb$tick())
  } else {
    opts <- NULL
  }
  parallel::clusterExport(cl = cl, varlist = c(
    "internal_selection", "logit_model_nonprobsvy", "start_fit", "get_method", "controlSel", "mle",
    "probit_model_nonprobsvy", "cloglog_model_nonprobsvy", "mm", "u_theta_beta_dr",
    "mu_hatDR"
  ))

  k <- 1:num_boot
  mu_hats <- foreach::`%dopar%`(
    obj = foreach::foreach(k = k, .combine = c, .options.snow = opts),
    ex = {
      strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
      # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)

      weights_strap <- weights[strap_nons]
      # weights_rand_strap <- weights_rand[strap_rand]

      # using svy package
      strap_rand_svy <- which(rep_weights[, k] != 0)
      weights_rand_strap_svy <- rep_weights[, k] * weights_rand
      # N_strap <- sum(weights_rand_strap_svy)
      # X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
      weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

      X_strap <- rbind(X_rand[strap_rand_svy, , drop = FALSE], X_nons[strap_nons, , drop = FALSE])
      y_strap <- c(y_rand[strap_rand_svy], y_nons[strap_nons])
      n_rand_strap <- nrow(X_rand[strap_rand_svy, , drop = FALSE])

      R_nons <- rep(1, n_nons)
      R_rand <- rep(0, n_rand_strap)
      R <- c(R_rand, R_nons)

      model_strap <- mm(
        X = X_strap,
        y = y_strap,
        weights = weights_strap,
        weights_rand = weights_strap_rand,
        R = R, # c(R[loc_nons][strap_nons], R[loc_rand][strap_rand]),
        n_nons = n_nons,
        n_rand = n_rand_strap,
        method_selection = method_selection,
        family = family_nonprobsvy,
        start_selection = start_selection,
        start_outcome = start_outcome,
        boot = TRUE
      )

      weights_nons_strap <- 1 / model_strap$selection$ps_nons
      N_nons <- sum(weights_strap * weights_nons_strap)
      N_rand <- sum(weights_strap_rand)

      mu_hatDR(
        y = y_nons[strap_nons],
        y_nons = model_strap$outcome$y_nons_pred,
        y_rand = model_strap$outcome$y_rand_pred,
        weights = weights_strap,
        weights_nons = weights_nons_strap,
        weights_rand = weights_strap_rand,
        N_nons = N_nons,
        N_rand = N_rand
      ) # DR estimator
    }
  )
  # mu_hat_boot <- mean(mu_hats)
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}
