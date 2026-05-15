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
                    case_weights,
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

  # One bootstrap iteration; returns a scalar y_mi_hat or stops on error.
  # Centralising the body lets single-core and multicore paths share logic
  # and ensures the same fix is applied to both.
  one_iter <- function(b) {
    if (!is.null(svydesign)) {
      strap_rand_svy <- which(rep_weights[, b] != 0)
      weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
      pop_size_strap <- sum(weights_rand_strap_svy)

      # Subset the replicate design through `[.svyrep.design` so that all
      # internal slots (variables, repweights, pweights, selfrep, degf) stay
      # consistent. Mutating slots by hand left `selfrep` at its original
      # length and caused `survey::svytotal()` to fail deterministically
      # with "(subscript) logical subscript too long" - which the previous
      # silent tryCatch around the bootstrap loop swallowed, causing the
      # loop counter never to advance and the call to hang for ever.
      svydesign_rep_b <- svydesign_rep[strap_rand_svy, ]

      strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons),
                               prob = 1 / case_weights)
      model_obj_b <- outcome_method(
        y_nons = y_nons[strap_nons],
        X_nons = X_nons[strap_nons, , drop = FALSE],
        X_rand = X_rand[strap_rand_svy, , drop = FALSE],
        svydesign = svydesign_rep_b,
        weights = case_weights[strap_nons],
        family_outcome = family_outcome,
        start_outcome = model_obj$coefficients,
        vars_selection = model_obj$vars_selection,
        pop_totals = pop_totals,
        pop_size = pop_size_strap,
        control_outcome = control_outcome,
        control_inference = control_inference,
        verbose = FALSE,
        se = FALSE
      )
    } else {
      strap_nons <- sample.int(replace = TRUE, n = NROW(X_nons),
                               prob = 1 / case_weights)
      model_obj_b <- outcome_method(
        y_nons = y_nons[strap_nons],
        X_nons = X_nons[strap_nons, , drop = FALSE],
        X_rand = X_rand,
        svydesign = svydesign,
        weights = case_weights[strap_nons],
        family_outcome = family_outcome,
        start_outcome = model_obj$coefficients,
        vars_selection = model_obj$vars_selection,
        pop_totals = pop_totals,
        pop_size = pop_size,
        control_outcome = control_outcome,
        control_inference = control_inference,
        verbose = FALSE,
        se = FALSE
      )
    }
    model_obj_b$y_mi_hat
  }

  # Single core processing
  if (control_inference$cores == 1) {

    if (verbose) {
      message("Single core bootstrap in progress...")
      pb_boot <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
    }

    # Allow at most one retry per replicate so a deterministic failure
    # surfaces as a proper error instead of looping forever. The previous
    # implementation used `while (b <= num_boot)` with a silent tryCatch,
    # which made deterministic failures invisible and could hang the
    # caller indefinitely.
    max_retries <- 1L
    for (b in seq_len(num_boot)) {
      attempt <- 0L
      repeat {
        res <- tryCatch(one_iter(b),
                        error = function(e) e)
        if (!inherits(res, "error")) {
          boot_obj[b] <- res
          break
        }
        attempt <- attempt + 1L
        if (verbose) {
          message(sprintf(
            "An error occurred in iteration %d (attempt %d/%d): %s",
            b, attempt, max_retries + 1L, conditionMessage(res)))
        }
        if (attempt > max_retries) {
          stop(sprintf(
            "Bootstrap iteration %d failed after %d attempts: %s",
            b, max_retries + 1L, conditionMessage(res)),
            call. = FALSE)
        }
      }
      if (verbose) utils::setTxtProgressBar(pb_boot, b)
    }
    if (verbose) close(pb_boot)
  } else {
    # Multicore processing
    if (verbose) message("Multicore bootstrap in progress...")

    cl <- parallel::makeCluster(control_inference$cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl = cl, varlist = NULL, envir = getNamespace("nonprobsvy"))

    boot_obj <- foreach::`%dopar%`(
      obj = foreach::foreach(b = 1:num_boot, .combine = c),
      ex = {
        one_iter(b)
      }
    )
  }


  # Return results
  return(boot_obj)
}
