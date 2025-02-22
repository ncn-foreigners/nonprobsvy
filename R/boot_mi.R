#' Bootstrap for the MI estimator
#' @param model_obj model object (from method_*)
#' @param y_nons target variable vector
#' @param X_nons auxiliary variable matrix (non-probability)
#' @param X_rand auxiliary variable matrix (probability sample)
#' @param weights case weights
#' @param pop_totals population totals
#' @param pop_size population size
#' @param family_outcome family for the method_* functions
#' @param control_outcome controls for the outcome
#' @param control_inference controls for the inference
#' @param verbose whether to print information
#' @noRd
boot_mi <- function(model_obj,
                    y_nons,
                    X_nons,
                    X_rand,
                    weights,
                    pop_totals,
                    pop_size,
                    family_outcome,
                    control_outcome,
                    control_inference,
                    verbose) {

  outcome_method <- switch(model_obj$model,
                           "glm" = method_glm,
                           "nn" = method_nn,
                           "pmm" = method_pmm,
                           "npar" = method_npar)

  # Initialize objects to store results
  num_boot <- control_inference$num_boot
  svydesign <- model_obj$svydesign

  boot_obj <- numeric(num_boot)

  # Prepare survey design replicates if available
  if (!is.null(svydesign)) {
    svydesign_rep <- survey::as.svrepdesign(svydesign,
                                            type = control_inference$rep_type,
                                            replicates = num_boot)
    rep_weights <- svydesign_rep$repweights$weights
  }

  # Single core processing
  if (control_inference$cores == 1) {

    if (verbose) {
      message("Single core bootstrap in progress...")
      pb_boot <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
    }

    b <- 1

    if (!is.null(svydesign)) {
      # Bootstrap for probability and non-probability samples
      while (b <= num_boot) {
        tryCatch(
          {
            # Probability part
            strap_rand_svy <- which(rep_weights[, b] != 0)
            weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
            pop_size_strap <- sum(weights_rand_strap_svy)

            # Workaround for as.svrepdesign
            svydesign_rep_b <- svydesign_rep
            svydesign_rep_b$variables <- svydesign_rep$variables[strap_rand_svy, ]
            svydesign_rep_b$repweights <- svydesign_rep$repweights[strap_rand_svy, ]
            svydesign_rep_b$pweights <- svydesign_rep$pweights[strap_rand_svy]

            # Non-probability part
            strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons), prob = 1 / weights)

            model_obj_b <- outcome_method(y_nons = y_nons[strap_nons],
                                          X_nons = X_nons[strap_nons, , drop = FALSE],
                                          X_rand = X_rand[strap_rand_svy, , drop = FALSE],
                                          svydesign = svydesign_rep_b,
                                          weights = weights[strap_nons],
                                          family_outcome = family_outcome,
                                          start_outcome = model_obj$coefficients,
                                          vars_selection = model_obj$vars_selection,
                                          pop_totals = pop_totals,
                                          pop_size = pop_size_strap,
                                          control_outcome = control_outcome,
                                          control_inference = control_inference,
                                          verbose = FALSE,
                                          se = FALSE)

            boot_obj[b] <- model_obj_b$y_mi_hat

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
        tryCatch(
          {
            # Non-probability part
            strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons), prob = 1 / weights)

            model_obj_b <- outcome_method(y_nons = y_nons[strap_nons],
                                          X_nons = X_nons[strap_nons, , drop = FALSE],
                                          X_rand = X_rand,
                                          svydesign = svydesign,
                                          weights = weights[strap_nons],
                                          family_outcome = family_outcome,
                                          start_outcome = model_obj$coefficients,
                                          vars_selection = model_obj$vars_selection,
                                          pop_totals = pop_totals,
                                          pop_size = pop_size,
                                          control_outcome = control_outcome,
                                          control_inference = control_inference,
                                          verbose = FALSE,
                                          se = FALSE)

            boot_obj[b] <- model_obj_b$y_mi_hat

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
        obj = foreach::foreach(b = 1:num_boot, .combine = c),
        ex = {
          strap_rand_svy <- which(rep_weights[, b] != 0)
          weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
          pop_size_strap <- sum(weights_rand_strap_svy)

          # Workaround for as.svrepdesign
          svydesign_rep_b <- svydesign_rep
          svydesign_rep_b$variables <- svydesign_rep$variables[strap_rand_svy, ]
          svydesign_rep_b$repweights <- svydesign_rep$repweights[strap_rand_svy, ]
          svydesign_rep_b$pweights <- svydesign_rep$pweights[strap_rand_svy]

          # Non-probability part
          strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons), prob = 1 / weights)

          model_obj_b <- outcome_method(y_nons = y_nons[strap_nons],
                                        X_nons = X_nons[strap_nons, , drop = FALSE],
                                        X_rand = X_rand[strap_rand_svy, , drop = FALSE],
                                        svydesign = svydesign_rep_b,
                                        weights = weights[strap_nons],
                                        family_outcome = family_outcome,
                                        start_outcome = model_obj$coefficients,
                                        vars_selection = model_obj$vars_selection,
                                        pop_totals = pop_totals,
                                        pop_size = pop_size_strap,
                                        control_outcome = control_outcome,
                                        control_inference = control_inference,
                                        verbose = FALSE,
                                        se = FALSE)
          model_obj_b$y_mi_hat
        }
      )
    } else {
      # Parallel bootstrap for non-probability samples only
      boot_obj <- foreach::`%dopar%`(
        obj = foreach::foreach(b = 1:num_boot, .combine = c),
        ex = {
          strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons), prob = 1 / weights)

          model_obj_b <- outcome_method(y_nons = y_nons[strap_nons],
                                        X_nons = X_nons[strap_nons, , drop = FALSE],
                                        X_rand = X_rand,
                                        svydesign = svydesign,
                                        weights = weights[strap_nons],
                                        family_outcome = family_outcome,
                                        start_outcome = model_obj$coefficients,
                                        vars_selection = model_obj$vars_selection,
                                        pop_totals = pop_totals,
                                        pop_size = pop_size,
                                        control_outcome = control_outcome,
                                        control_inference = control_inference,
                                        verbose = FALSE,
                                        se = FALSE)

          model_obj_b$y_mi_hat
        }
      )
    }
  }


  # Return results
  return(boot_obj)
}
