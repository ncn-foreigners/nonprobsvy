## bootstrap for the IPW estimator
boot_ipw <- function(X_rand,
                     X_nons,
                     svydesign,
                     case_weights,
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
                     gee_h_fun,
                     rep_type,
                     maxit,
                     control_inference,
                     control_selection,
                     verbose,
                     pop_size,
                     pop_totals,
                     ...) {

  if (!is.null(weights_rand)) N <- sum(weights_rand)

  estimation_method <- est_method_ipw(est_method)

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

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
    rep_weights <- survey::as.svrepdesign(design = svydesign,
                                          type = rep_type,
                                          replicates = num_boot)$repweights$weights # TODO customise to calibrated svydesign
    while (k <= num_boot) {
      tryCatch(
        {
          strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / case_weights)

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

          est_method_obj <- estimation_method$estimation_model(
            model = estimation_method$model_selection(
              X = X,
              X_nons = X_nons_strap,
              X_rand = X_rand_strap,
              weights = case_weights[strap_nons],
              weights_rand = weights_strap_rand,
              R = R,
              method_selection = method_selection,
              optim_method = optim_method,
              gee_h_fun = gee_h_fun,
              est_method = est_method,
              maxit = maxit,
              start = start_selection,
              control_selection = control_selection,
              verbose = verbose,
              varcov = TRUE
            ),
            method_selection = method_selection
          )

          ps_nons <- est_method_obj$ps_nons
          weights_nons <- 1 / ps_nons
          N_est_nons <- ifelse(is.null(pop_size), sum(case_weights[strap_nons] * weights_nons), pop_size)

          for (l in 1:mu_len) {
            mu_hats_boot[k, l] <- mu_hatIPW(
              y = ys[[l]][strap_nons],
              weights = case_weights[strap_nons],
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
            message(info)
          }
        }
      )
    }
  } else {
    while (k <= num_boot) {
      tryCatch(
        {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / case_weights)

          X_strap <- X_nons[strap, , drop = FALSE]
          R_strap <- R[strap]
          weights_strap <- case_weights[strap]

          h_object_strap <- theta_h_estimation(
            R = R_strap,
            X = X_strap,
            weights_rand = NULL,
            weights = weights_strap,
            gee_h_fun = gee_h_fun,
            method_selection = method_selection,
            maxit = maxit,
            pop_totals = pop_totals,
            nleqslv_method = control_selection$nleqslv_method,
            nleqslv_global = control_selection$nleqslv_global,
            nleqslv_xscalm = control_selection$nleqslv_xscalm
          )
          theta_hat_strap <- h_object_strap$theta_h
          ps_nons <- inv_link(tcrossprod(theta_hat_strap,X_strap))

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
            message(info)
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

# Multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
boot_ipw_multicore <- function(X_rand,
                               X_nons,
                               svydesign,
                               case_weights,
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
                               gee_h_fun,
                               maxit,
                               control_selection,
                               control_inference,
                               cores,
                               pop_size,
                               pop_totals,
                               verbose,
                               ...) {

  if (!is.null(weights_rand)) N <- sum(weights_rand)

  estimation_method <- est_method_ipw(est_method)

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  inv_link <- method$make_link_inv
  rep_type <- control_inference$rep_type

  mu_len <- length(mu_hats)
  mu_hats_boot <- numeric(length = num_boot * mu_len)
  boot_vars <- numeric(length = mu_len)

  if (verbose) {
    message("Multicore bootstrap in progress....")
  }

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl = cl,
                          varlist = c("method_ps", "est_method_ipw", "control_sel",
                                      "mu_hatIPW", "theta_h_estimation"),
                          envir = getNamespace("nonprobsvy"))

  rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights

  k <- 1:num_boot
  mu_hats_boot <- foreach::`%dopar%`(
    obj = foreach::foreach(k = k, .combine = c),
    ex = {
      if (is.null(pop_totals)) {
        strap_nons <- sample.int(replace = TRUE, n = n_nons, prob = 1 / case_weights)

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

        est_method_obj <- estimation_method$estimation_model(
          model = estimation_method$model_selection(
            X = X,
            X_nons = X_nons_strap,
            X_rand = X_rand_strap,
            weights = case_weights[strap_nons],
            weights_rand = weights_strap_rand,
            R = R,
            method_selection = method_selection,
            optim_method = optim_method,
            gee_h_fun = gee_h_fun,
            est_method = est_method,
            maxit = maxit,
            start = start_selection,
            control_selection = control_selection,
            verbose = verbose,
            varcov = TRUE
          ),
          method_selection = method_selection
        )

        ps_nons <- est_method_obj$ps_nons
        weights_nons <- 1 / ps_nons
        N_est_nons <- ifelse(is.null(pop_size), sum(case_weights[strap_nons] * weights_nons), pop_size)
        mu_hats_this_boot <- numeric(mu_len)

        for (l in 1:mu_len) {
          mu_hats_this_boot[l] <- mu_hatIPW(
            y = ys[[l]][strap_nons],
            weights = case_weights[strap_nons],
            weights_nons = weights_nons,
            N = N_est_nons
          ) # IPW estimator
        }
        mu_hats_this_boot
      } else {
        strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / case_weights)
        X_strap <- X_nons[strap, , drop = FALSE]
        R_strap <- R[strap]
        weights_strap <- case_weights[strap]

        h_object_strap <- theta_h_estimation(
          R = R_strap,
          X = X_strap,
          weights_rand = NULL,
          weights = weights_strap,
          gee_h_fun = gee_h_fun,
          method_selection = method_selection,
          maxit = maxit,
          pop_totals = pop_totals,
          start = start_selection,
          nleqslv_method = control_selection$nleqslv_method,
          nleqslv_global = control_selection$nleqslv_global,
          nleqslv_xscalm = control_selection$nleqslv_xscalm
        )
        theta_hat_strap <- h_object_strap$theta_h
        ps_nons <- inv_link(tcrossprod(theta_hat_strap,X_strap))

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
        mu_hats_boot
      }
    }
  )
  mu_hats_boot <- matrix(mu_hats_boot, nrow = num_boot, ncol = mu_len, byrow = TRUE)
  for (l in 1:mu_len) {
    boot_vars[l] <- 1 / (num_boot - 1) * sum((mu_hats_boot[, l] - mu_hats[l])^2)
  }
  list(
    var = boot_vars,
    # mu = mu_hats_boot_means,
    stat = mu_hats_boot
  )
}
