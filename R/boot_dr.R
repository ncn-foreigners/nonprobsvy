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

#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
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
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    ## progress bar
    # if (verbose) {
    #   pb <- progress::progress_bar$new(total = num_boot)
    #   opts <- list(progress = \(n) pb$tick())
    # } else {
    #   opts <- NULL
    # }
    parallel::clusterExport(cl = cl, varlist = c(
      "internal_selection", "internal_outcome", "logit_model_nonprobsvy", "start_fit", "get_method", "controlSel", "theta_h_estimation",
      "mle", "mu_hatDR", "probit_model_nonprobsvy", "cloglog_model_nonprobsvy", "glm_nonprobsvy", "nn_nonprobsvy", "pmm_nonprobsvy",
      "gaussian_nonprobsvy", "poisson_nonprobsvy", "binomial_nonprobsvy", "nonprobMI_fit", "controlOut"
    ))
    if (is.null(pop_totals)) {
      N <- sum(weights_rand)

      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c),
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
