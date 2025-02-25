#' Bootstrap for the DR estimator
#' @noRd
boot_dr <- function(selection,
                    outcome,
                    target,
                    data,
                    svydesign,
                    pop_totals,
                    pop_means,
                    pop_size,
                    method_selection,
                    method_outcome,
                    family_outcome,
                    subset,
                    strata,
                    weights,
                    na_action,
                    control_selection,
                    control_outcome,
                    control_inference,
                    start_outcome,
                    start_selection,
                    verbose,
                    pop_size_fixed) {

  # Initialize objects to store results
  num_boot <- control_inference$num_boot

  ## add bootstrap with variables combination

  if (!is.null(svydesign)) {
    svydesign_rep <- survey::as.svrepdesign(svydesign,
                                            type = control_inference$rep_type,
                                            replicates = num_boot)
    rep_weights <- svydesign_rep$repweights$weights
  }

  # Single core processing
  if (control_inference$cores == 1) {

    boot_obj <- matrix(0, ncol = length(all.vars(target)), nrow = num_boot)

    if (verbose) {
      message("Single core bootstrap in progress...")
      pb_boot <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
    }

    b <- 1

    if (!is.null(svydesign)) {


      # Bootstrap for probability and non-probability samples
      while (b <= num_boot) {

        ## probability part
        strap_rand_svy <- which(rep_weights[, b] != 0)
        weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
        pop_size_strap <- sum(weights_rand_strap_svy)

        data_prob <- svydesign$variables[strap_rand_svy, ]
        data_prob$weight <- weights_rand_strap_svy[strap_rand_svy]

        svyd_call <- as.list(svydesign$call)
        svyd_call[[1]] <- NULL
        svyd_call$ids <- as.formula(svyd_call$ids)
        svyd_call$weights <- as.formula(svyd_call$weights)
        svyd_call$strata <- as.formula(svyd_call$strata)
        svyd_call$data <- as.name("data_prob")

        # Method 1: Using do.call
        svydesign_b <- do.call(survey::svydesign, svyd_call)

        strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / weights)


        tryCatch(
          {
            results_ipw_b <- nonprob_ipw(selection = selection,
                                       target = target,
                                       data = data[strap_nons, ],
                                       svydesign = svydesign_b,
                                       pop_totals = NULL,
                                       pop_means = NULL,
                                       pop_size = NULL,
                                       method_selection = method_selection,
                                       subset = subset,
                                       strata = strata,
                                       weights = weights[strap_nons],
                                       na_action = na_action,
                                       control_selection = control_selection,
                                       control_inference = control_inference,
                                       start_selection = start_selection,
                                       verbose = FALSE,
                                       se = FALSE,
                                       pop_size_fixed=pop_size_fixed)
            ## estimate the mi
            results_mi_b <- nonprob_mi(outcome = outcome,
                                       data = data[strap_nons, ],
                                     svydesign = svydesign_b,
                                     pop_totals = NULL,
                                     pop_means = NULL,
                                     pop_size = NULL,
                                     method_outcome = method_outcome,
                                     family_outcome = family_outcome,
                                     subset = subset,
                                     strata = strata,
                                     weights = weights[strap_nons],
                                     na_action = na_action,
                                     control_outcome = control_outcome,
                                     control_inference = control_inference,
                                     start_outcome = start_outcome,
                                     verbose = FALSE,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)

            boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                      y_resid = do.call("cbind", results_mi_b$ys_resid),
                                      weights = weights[strap_nons],
                                      weights_nons = results_ipw_b$ipw_weights,
                                      N_nons = sum(results_ipw_b$ipw_weights))

            if (verbose) {
              utils::setTxtProgressBar(pb_boot, b)
            }
            b <- b + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", b, " iteration: ", e$message, sep = "")
              message(info)
            }
          }
        )
      }
    } else {
      # Bootstrap for non-probability samples only
      while (b <= num_boot) {

        strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / weights)

        tryCatch(
          {
            # Non-probability part
            results_ipw_b <- nonprob_ipw(selection = selection,
                                       target = target,
                                       data = data[strap_nons, ],
                                       svydesign = svydesign,
                                       pop_totals = pop_totals,
                                       pop_means = pop_means,
                                       pop_size = pop_size,
                                       method_selection = method_selection,
                                       subset = subset,
                                       strata = strata,
                                       weights = weights[strap_nons],
                                       na_action = na_action,
                                       control_selection = control_selection,
                                       control_inference = control_inference,
                                       start_selection = start_selection,
                                       verbose = FALSE,
                                       se = FALSE,
                                       pop_size_fixed=pop_size_fixed)
            ## estimate the mi
            results_mi_b <- nonprob_mi(outcome = outcome,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign,
                                     pop_totals = pop_totals,
                                     pop_means = pop_means,
                                     pop_size = pop_size,
                                     method_outcome = method_outcome,
                                     family_outcome = family_outcome,
                                     subset = subset,
                                     strata = strata,
                                     weights = weights[strap_nons],
                                     na_action = na_action,
                                     control_outcome = control_outcome,
                                     control_inference = control_inference,
                                     start_outcome = start_outcome,
                                     verbose = FALSE,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)

            boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                      y_resid = do.call("cbind", results_mi_b$ys_resid),
                                      weights = weights[strap_nons],
                                      weights_nons = results_ipw_b$ipw_weights,
                                      N_nons = sum(results_ipw_b$ipw_weights))

            if (verbose) {
              utils::setTxtProgressBar(pb_boot, b)
            }
            b <- b + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", b, " iteration: ", e$message, sep = "")
              message(info)
            }
          }
        )
      }
    }
  } else {
    # Multicore processing
    if (verbose) message("Multicore bootstrap in progress...")

    cl <- parallel::makeCluster(control_inference$cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl = cl, varlist = NULL, envir = getNamespace("nonprobsvy"))

    if (!is.null(svydesign)) {
      # Parallel bootstrap for probability and non-probability samples
      boot_obj <- foreach::`%dopar%`(
        obj = foreach::foreach(b = 1:num_boot, .combine = rbind),
        ex = {

          strap_rand_svy <- which(rep_weights[, b] != 0)
          weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
          pop_size_strap <- sum(weights_rand_strap_svy)

          data_prob <- svydesign$variables[strap_rand_svy, ]
          data_prob$weight <- weights_rand_strap_svy[strap_rand_svy]

          svyd_call <- as.list(svydesign$call)
          svyd_call[[1]] <- NULL
          svyd_call$ids <- as.formula(svyd_call$ids)
          svyd_call$weights <- as.formula(svyd_call$weights)
          svyd_call$strata <- as.formula(svyd_call$strata)
          svyd_call$data <- as.name("data_prob")

          # Method 1: Using do.call
          svydesign_b <- do.call(survey::svydesign, svyd_call)

          strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / weights)

          results_ipw_b <- nonprob_ipw(selection = selection,
                                     target = target,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign_b,
                                     pop_totals = NULL,
                                     pop_means = NULL,
                                     pop_size = NULL,
                                     method_selection = method_selection,
                                     subset = subset,
                                     strata = strata,
                                     weights = weights[strap_nons],
                                     na_action = na_action,
                                     control_selection = control_selection,
                                     control_inference = control_inference,
                                     start_selection = start_selection,
                                     verbose = verbose,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)
          ## estimate the mi
          results_mi_b <- nonprob_mi(outcome = outcome,
                                     data = data[strap_nons, ],
                                   svydesign = svydesign_b,
                                   pop_totals = NULL,
                                   pop_means = NULL,
                                   pop_size = NULL,
                                   method_outcome = method_outcome,
                                   family_outcome = family_outcome,
                                   subset = subset,
                                   strata = strata,
                                   weights = weights[strap_nons],
                                   na_action = na_action,
                                   control_outcome = control_outcome,
                                   control_inference = control_inference,
                                   start_outcome = start_outcome,
                                   verbose = verbose,
                                   se = FALSE,
                                   pop_size_fixed=pop_size_fixed)

          boot_obj_b <- mu_hatDR(y_hat = results_mi_b$output$mean,
                   y_resid = do.call("cbind", results_mi_b$ys_resid),
                   weights = weights[strap_nons],
                   weights_nons = results_ipw_b$ipw_weights,
                   N_nons = sum(results_ipw_b$ipw_weights))

          as.matrix(boot_obj_b)
        }
      )
    } else {
      # Parallel bootstrap for non-probability samples only
      boot_obj <- foreach::`%dopar%`(
        obj = foreach::foreach(b = 1:num_boot, .combine = rbind),
        ex = {

          strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / weights)

          results_ipw_b <- nonprob_ipw(selection = selection,
                                     target = target,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign,
                                     pop_totals = pop_totals,
                                     pop_means = pop_means,
                                     pop_size = pop_size,
                                     method_selection = method_selection,
                                     subset = subset,
                                     strata = strata,
                                     weights = weights[strap_nons],
                                     na_action = na_action,
                                     control_selection = control_selection,
                                     control_inference = control_inference,
                                     start_selection = start_selection,
                                     verbose = verbose,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)
          ## estimate the mi
          results_mi_b <- nonprob_mi(outcome = outcome,
                                     data = data[strap_nons, ],
                                   svydesign = svydesign,
                                   pop_totals = pop_totals,
                                   pop_means = pop_means,
                                   pop_size = pop_size,
                                   method_outcome = method_outcome,
                                   family_outcome = family_outcome,
                                   subset = subset,
                                   strata = strata,
                                   weights = weights[strap_nons],
                                   na_action = na_action,
                                   control_outcome = control_outcome,
                                   control_inference = control_inference,
                                   start_outcome = start_outcome,
                                   verbose = verbose,
                                   se = FALSE,
                                   pop_size_fixed=pop_size_fixed)

          boot_obj_b <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                 y_resid = do.call("cbind", results_mi_b$ys_resid),
                                 weights = weights[strap_nons],
                                 weights_nons = results_ipw_b$ipw_weights,
                                 N_nons = sum(results_ipw_b$ipw_weights))
          as.matrix(boot_obj_b)
        }
      )
    }
  }


  # Return results
  return(boot_obj)
}
