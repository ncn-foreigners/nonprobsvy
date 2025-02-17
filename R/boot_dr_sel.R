boot_dr_sel <- function(X,
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
#' @importFrom doParallel registerDoParallel
boot_dr_sel_multicore <- function(X,
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
                                  verbose) {
  mu_hats <- vector(mode = "numeric", length = num_boot)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- X[loc_nons, , drop = FALSE]
  X_rand <- X[loc_rand, , drop = FALSE]
  y_nons <- y[loc_nons]
  y_rand <- y[loc_rand]

  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

  if (verbose) message("Multicores bootstrap in progress..")

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl = cl, varlist = c(
    "internal_selection", "model_ps", "start_fit", "get_method", "control_sel", "mle", "mm", "u_theta_beta_dr",
    "mu_hatDR"
  ))

  k <- 1:num_boot
  mu_hats <- foreach::`%dopar%`(
    obj = foreach::foreach(k = k, .combine = c),
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
