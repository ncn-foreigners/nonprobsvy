#' @useDynLib nonprobsvy
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats weighted.mean
#' @importFrom RANN nn2
#' @importFrom ncvreg cv.ncvreg
#' @importFrom stats terms
#' @importFrom stats reformulate
#' @importFrom survey svytotal
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @rdname nonprob
#' @export
nonprob_mi <- function(outcome,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_outcome,
                       family_outcome = "gaussian",
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_outcome,
                       control_inference,
                       start_outcome,
                       verbose,
                       x,
                       y,
                       se,
                       ...) {
  var_selection <- control_inference$vars_selection
  outcome_init <- outcome
  outcomes <- make_outcomes(outcome)
  output <- list()
  ys <- list()
  outcome_list <- list()
  if (se) {
    confidence_interval <- list()
    SE_values <- list()
  } else {
    confidence_interval <- NULL
    SE_values <- NULL
  }
  num_boot <- control_inference$num_boot
  if (method_outcome == "pmm" & (!is.null(pop_totals) | !is.null(pop_means))) {
    control_inference$var_method <- "bootstrap"
    message("Bootstrap variance only, analytical version during implementation.")
  }
  if (control_inference$var_method == "bootstrap") {
    stat <- matrix(nrow = control_inference$num_boot, ncol = outcomes$l)
  }
  terms_obj <- terms(outcome)
  if (is.null(pop_totals) && !is.null(svydesign)) {
    prob_totals <- svytotal(reformulate(attr(terms_obj, "term.labels")), svydesign)
    prob_pop_totals <- c(sum(weights(svydesign)), prob_totals)
  }

  for (k in 1:outcomes$l) {
    if (control_outcome$pmm_k_choice == "min_var" & method_outcome == "pmm") {
      # This can be programmed a lot better possibly with custom method outcome that would
      # store previous k-pmm model and omit the last estimation

      ## TODO:: right now this only goes forward not backwards
      var_prev <- Inf
      cond <- TRUE
      kk <- control_outcome$k - 1
      while (cond) {
        outcome_model_data <- make_model_frame(formula = outcome, data = data, svydesign = svydesign)
        X_nons <- outcome_model_data$X_nons
        X_rand <- outcome_model_data$X_rand
        nons_names <- outcome_model_data$nons_names
        y_nons <- outcome_model_data$y_nons

        R_nons <- rep(1, nrow(X_nons))
        R_rand <- rep(0, nrow(X_rand))
        R <- c(R_nons, R_rand)

        loc_nons <- which(R == 1)
        loc_rand <- which(R == 0)

        n_nons <- nrow(X_nons)
        n_rand <- nrow(X_rand)
        X <- rbind(X_nons, X_rand)

        ps_rand <- svydesign$prob
        weights_rand <- 1 / ps_rand
        N_est_rand <- sum(weights_rand)

        kk <- kk + 1
        method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
        ## estimation

        MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
        model_obj <- MethodOutcome(
          outcome = outcome,
          data = data,
          weights = weights,
          family_outcome = family_outcome,
          start_outcome = start_outcome,
          X_nons = X_nons,
          y_nons = y_nons,
          X_rand = X_rand,
          control = control_outcome,
          n_nons = n_nons,
          n_rand = n_rand,
          model_frame = outcome_model_data$model_frame_rand,
          vars_selection = control_inference$vars_selection,
          pop_totals = pop_totals
        )
        y_rand_pred <- model_obj$y_rand_pred
        y_nons_pred <- model_obj$y_nons_pred
        # parameters <- model_obj$parameters
        outcome_list[[k]] <- model_obj$model

        # updating probability sample by adding y_hat variable
        svydesign1 <- stats::update(svydesign,
          y_hat_MI = y_rand_pred
        )
        mu_hat <- weighted.mean(y_rand_pred, w = weights_rand)

        var_obj <- internal_varMI(
          svydesign = svydesign1,
          X_nons = X_nons,
          X_rand = X_rand,
          y = y_nons,
          y_pred = y_nons_pred,
          weights_rand = weights_rand,
          method = method_outcome,
          n_rand = n_rand,
          n_nons = n_nons,
          N = N_est_rand,
          family = family_outcome,
          model_obj = model_obj,
          pop_totals = pop_totals,
          k = control_outcome$k,
          pmm_match_type = control_outcome$pmm_match_type,
          nn_exact_se = control_inference$nn_exact_se,
          pmm_reg_engine = control_outcome$pmm_reg_engine,
          pi_ij = control_inference$pi_ij
        )

        var_nonprob <- var_obj$var_nonprob
        var_prob <- var_obj$var_prob

        se_nonprob <- sqrt(var_nonprob)
        se_prob <- sqrt(var_prob)
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
        # variance
        var_now <- var_nonprob + var_prob
        cond <- var_prev > var_now
        var_prev <- var_now
      }
      control_outcome$k <- kk - 1
      svydesign1 <- NULL # freeing up memmory
    }

    if (is.null(pop_totals) && !is.null(svydesign)) {
      pop_totals_sel <- pop_totals
      outcome <- outcomes$outcome[[k]]

      # model for outcome formula
      outcome_model_data <- make_model_frame(formula = outcome, data = data, svydesign = svydesign)
      X_nons <- outcome_model_data$X_nons
      X_rand <- outcome_model_data$X_rand
      nons_names <- outcome_model_data$nons_names
      y_nons <- outcome_model_data$y_nons

      R_nons <- rep(1, nrow(X_nons))
      R_rand <- rep(0, nrow(X_rand))
      R <- c(R_nons, R_rand)

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(X_nons)
      n_rand <- nrow(X_rand)
      X <- rbind(X_nons, X_rand)

      ps_rand <- svydesign$prob
      weights_rand <- 1 / ps_rand
      N_est_rand <- sum(weights_rand)

      ########### WORKING VERSION

      if (var_selection == TRUE) {
        # TODO add variables randomization
        nlambda <- control_outcome$nlambda
        beta <- ncvreg::cv.ncvreg(
          X = X_nons[, -1, drop = FALSE],
          y = y_nons,
          penalty = control_outcome$penalty,
          family = family_outcome,
          trace = verbose,
          nfolds = control_outcome$nfolds,
          nlambda = nlambda,
          gamma = switch(control_outcome$penalty,
            SCAD = control_outcome$a_SCAD,
            control_outcome$a_MCP
          ),
          lambda_min = control_outcome$lambda_min,
          eps = control_outcome$epsilon
        )

        beta_est <- beta$fit$beta[, beta$min]
        beta_selected <- which(abs(beta_est) != 0) - 1
        beta_est <- beta_est[beta$fit$beta[, beta$min] != 0]
        cve_outcome <- beta$cve
        lambda_outcome <- beta$lambda
        lambda_min_outcome <- beta$lambda.min

        X_design <- as.matrix(X[, beta_selected + 1, drop = FALSE])
        # colnames(X_design) <- c("(Intercept)", colnames(Xsel))
        X_rand <- X_design[loc_rand, , drop = FALSE]
        X_nons <- X_design[loc_nons, , drop = FALSE]
      }
      ################

      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      ## estimation

      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(
        outcome = outcome,
        data = data,
        weights = weights,
        family_outcome = family_outcome,
        start_outcome = start_outcome,
        X_nons = X_nons,
        y_nons = y_nons,
        X_rand = X_rand,
        control = control_outcome,
        n_nons = n_nons,
        n_rand = n_rand,
        model_frame = outcome_model_data$model_frame_rand,
        vars_selection = control_inference$vars_selection,
        pop_totals = pop_totals
      )
      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      # parameters <- model_obj$parameters
      outcome_list[[k]] <- model_obj$model
      outcome_list[[k]]$model_frame <- outcome_model_data$model_frame_rand

      # updating probability sample by adding y_hat variable
      svydesign <- stats::update(svydesign,
        y_hat_MI = y_rand_pred
      )
      mu_hat <- weighted.mean(y_rand_pred, w = weights_rand)
    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {
      if (!is.null(pop_totals)) pop_size <- N_est_rand <- pop_totals[1]
      if (!is.null(pop_means)) { # TO consider
        if (!is.null(pop_size)) {
          pop_totals <- c(pop_size, pop_size * pop_means)
          names(pop_totals) <- c("(Intercept)", names(pop_means))
        } else {
          stop("The `pop_size` argument must be specified when the `pop_means` argument is provided.")
        }
      }
      Model <- make_model_frame(formula = outcome, data = data, pop_totals = pop_totals)

      X_nons <- Model$X_nons
      X_rand <- Model$X_rand
      y_nons <- Model$y_nons
      n_rand <- 0
      weights_rand <- NULL
      n_nons <- nrow(X_nons)
      R <- rep(1, n_nons)

      ############## WORKING VERSION
      if (var_selection == TRUE) {
        nlambda <- control_outcome$nlambda
        beta <- ncvreg::cv.ncvreg(
          X = Model$X_nons[, -1, drop = FALSE],
          y = y_nons,
          penalty = control_outcome$penalty,
          family = family_outcome,
          trace = verbose,
          nfolds = control_outcome$nfolds,
          nlambda = nlambda,
          gamma = switch(control_outcome$penalty,
            SCAD = control_outcome$a_SCAD,
            control_outcome$a_MCP
          ),
          lambda_min = control_outcome$lambda_min,
          eps = control_outcome$epsilon
        )

        beta_est <- beta$fit$beta[, beta$min]
        beta_selected <- which(abs(beta_est) != 0) - 1
        beta_est <- beta_est[beta$fit$beta[, beta$min] != 0]
        cve_outcome <- beta$cve
        lambda_outcome <- beta$lambda
        lambda_min_outcome <- beta$lambda.min

        X_nons <- Model$X_nons[, beta_selected + 1, drop = FALSE]
        X <- X_nons
        pop_totals <- pop_totals[beta_selected + 1]
      }
      prob_pop_totals <- pop_totals
      #######################

      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(
        outcome = outcome,
        data = data,
        weights = weights,
        family_outcome = family_outcome,
        start_outcome = start_outcome,
        X_nons = X_nons,
        y_nons = y_nons,
        X_rand = X_rand,
        control = control_outcome,
        n_nons = n_nons,
        n_rand = n_rand,
        model_frame = Model$model_frame_rand,
        vars_selection = control_inference$vars_selection,
        pop_totals = pop_totals
      )

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      parameters <- model_obj$parameters
      outcome_list[[k]] <- model_obj$model
      outcome_list[[k]]$model_frame <- Model$model_frame_rand
      mu_hat <- y_rand_pred
    } else {
      stop("Specify one of the `svydesign`, `pop_totals' or `pop_means' arguments, not all.")
    }

    if (isTRUE(attr(model_obj$model, "method") == "pmm") & !(control_inference$nn_exact_se)) {
      # if not nn_exact_se then this can be dropped
      model_obj$model$glm_obj <- NULL
    }

    ys[[k]] <- as.numeric(y_nons)
    if (se) {
      # design based variance estimation based on approximations of the second-order inclusion probabilities
      if (control_inference$var_method == "analytic") { # consider move variance implementation to internals
        var_obj <- internal_varMI(
          svydesign = svydesign,
          X_nons = X_nons,
          X_rand = X_rand,
          y = y_nons,
          y_pred = y_nons_pred,
          weights_rand = weights_rand,
          method = method_outcome,
          n_rand = n_rand,
          n_nons = n_nons,
          N = N_est_rand,
          family = family_outcome,
          model_obj = model_obj,
          pop_totals = pop_totals,
          # we should probably just pass full control list
          k = control_outcome$k,
          pmm_match_type = control_outcome$pmm_match_type,
          nn_exact_se = control_inference$nn_exact_se,
          pmm_reg_engine = control_outcome$pmm_reg_engine,
          pi_ij = control_inference$pi_ij
        )

        var_nonprob <- var_obj$var_nonprob
        var_prob <- var_obj$var_prob

        se_nonprob <- sqrt(var_nonprob)
        se_prob <- sqrt(var_prob)
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
        # variance
        var <- var_nonprob + var_prob
      } else if (control_inference$var_method == "bootstrap") { # TODO for pop_totals
        # bootstrap variance
        if (control_inference$cores > 1) {
          boot_obj <- boot_mi_multicore(
            X_rand = X_rand,
            X_nons = X_nons,
            weights = weights,
            y = y_nons,
            family_outcome = family_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights_rand,
            mu_hat,
            svydesign,
            model_obj = model_obj,
            rep_type = control_inference$rep_type,
            method = method_outcome,
            control_outcome = control_outcome,
            control_inference = control_inference,
            pop_totals = pop_totals,
            cores = control_inference$cores,
            verbose = verbose
          )
        } else {
          boot_obj <- boot_mi(
            X_rand = X_rand,
            X_nons = X_nons,
            weights = weights,
            y = y_nons,
            family_outcome = family_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights_rand,
            mu_hat,
            svydesign,
            model_obj = model_obj,
            rep_type = control_inference$rep_type,
            method = method_outcome,
            control_outcome = control_outcome,
            control_inference = control_inference,
            pop_totals = pop_totals,
            verbose = verbose
          )
        }
        var <- boot_obj$var
        stat[, k] <- boot_obj$stat
        comp3_stat <- boot_obj$comp3_stat
        # mu_hat <- boot_obj$mu
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(
          nonprob = NA,
          prob = NA
        ))))
      } else {
        stop("Invalid `var_method` for the variance estimation.")
      }
      SE <- sqrt(var)
      alpha <- control_inference$alpha
      z <- stats::qnorm(1 - alpha / 2)
      # confidence interval based on the normal approximation
      confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(
        lower_bound = mu_hat - z * SE,
        upper_bound = mu_hat + z * SE
      ))))
    } else {
      SE <- NA
      confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(
        lower_bound = NA,
        upper_bound = NA
      ))))
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = NA, prob = NA))))
    }

    X <- rbind(X_rand, X_nons) # joint model matrix
    # if (is.null(pop_size))
    pop_size <- N_est_rand # estimated pop_size

    output[[k]] <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = SE))))
    outcome_list[[k]]$method <- method_outcome
    if (control_inference$vars_selection == TRUE) {
      outcome_list[[k]]$cve <- cve_outcome
    } else {
      NULL
    }
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f
  names(outcome_list) <- outcomes$f
  names(pop_size) <- "pop_size"
  names(ys) <- all.vars(outcome_init[[2]])
  names(prob_pop_totals) <- colnames(X_nons)

  boot_sample <- if (control_inference$var_method == "bootstrap" & control_inference$keep_boot) {
    list(stat = stat, comp2 = boot_obj$comp2)
  } else {
    NULL
  }
  if (!is.null(boot_sample) & is.matrix(boot_sample)) colnames(boot_sample) <- names(ys)

  structure(
    list(
      data = data,
      X = if (isTRUE(x)) X else NULL,
      y = if (isTRUE(y)) ys else NULL,
      R = R,
      prob = NULL,
      weights = NULL,
      control = list(
        control_outcome = control_outcome,
        control_inference = control_inference
      ),
      output = output,
      SE = SE_values,
      confidence_interval = confidence_interval,
      nonprob_size = n_nons,
      prob_size = n_rand,
      pop_size = pop_size,
      pop_totals = prob_pop_totals,
      outcome = outcome_list,
      selection = NULL,
      boot_sample = boot_sample,
      svydesign = if (is.null(pop_totals)) svydesign else NULL # TODO to customize if pop_totals only
    ),
    class = c("nonprobsvy", "nonprobsvy_mi")
  )
}
