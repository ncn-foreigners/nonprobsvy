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
#' @import Rcpp
#' @importFrom Rcpp evalCpp

nonprobMI <- function(outcome,
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
                      control_inference = controlInf(var_method = "analytic"),
                      start_outcome,
                      verbose,
                      x,
                      y,
                      se,
                      ...) {
  var_selection <- control_inference$vars_selection
  outcome_init <- outcome
  outcomes <- ff(outcome)
  output <- list()
  ys <- list()
  OutcomeList <- list()
  if (se) {
    confidence_interval <- list()
    SE_values <- list()
  } else {
    confidence_interval <- NULL
    SE_values <- NULL
  }
  num_boot <- control_inference$num_boot
  if (method_outcome == "pmm") {
    # control_inference$var_method <- "bootstrap"
    # message("Bootstrap variance only, analytical version during implementation.")
  }
  if (control_inference$var_method == "bootstrap") {
    stat <- matrix(nrow = control_inference$num_boot, ncol = outcomes$l)
  }
  for (k in 1:outcomes$l) {
    if (is.null(pop_totals) && !is.null(svydesign)) {
      pop_totals_sel <- pop_totals
      outcome <- outcomes$outcome[[k]]

      # model for outcome formula
      OutcomeModel <- model_frame(formula = outcome, data = data, svydesign = svydesign)
      X_nons <- OutcomeModel$X_nons
      X_rand <- OutcomeModel$X_rand
      nons_names <- OutcomeModel$nons_names
      y_nons <- OutcomeModel$y_nons

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
        X_rand <- X_design[loc_rand, ]
        X_nons <- X_design[loc_nons, ]
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
        model_frame = OutcomeModel$model_frame_rand,
        vars_selection = control_inference$vars_selection,
        pop_totals = pop_totals
      )
      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      # parameters <- model_obj$parameters
      OutcomeList[[k]] <- model_obj$model

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
          stop("pop_size must be defined when estimating with pop_means.")
        }
      }
      Model <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)

      X_nons <- Model$X_nons
      X_rand <- Model$X_rand
      y_nons <- Model$y_nons
      n_rand <- 0
      weights_rand <- NULL
      n_nons <- nrow(X_nons)

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
        model_frame = OutcomeModel$model_frame_rand,
        vars_selection = control_inference$vars_selection,
        pop_totals = pop_totals
      )

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      parameters <- model_obj$parameters
      OutcomeList[[k]] <- model_obj$model
      mu_hat <- y_rand_pred
    } else {
      stop("Please, provide svydesign object or pop_totals/pop_means.")
    }

    if (isTRUE(attr(model_obj$model, "method") == "pmm") & !(control_inference$pmm_exact_se)) {
      # if not pmm_exact_se then this can be dropped
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
          k = control_outcome$k,
          predictive_match = control_outcome$predictive_match,
          pmm_exact_se = control_inference$pmm_exact_se
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
          boot_obj <- bootMI_multicore(X_rand,
            X_nons,
            weights,
            y_nons,
            family_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights_rand,
            mu_hat,
            svydesign,
            rep_type = control_inference$rep_type,
            method = method_outcome,
            control_outcome = control_outcome,
            control_inference = control_inference,
            pop_totals = pop_totals,
            cores = control_inference$cores,
            verbose = verbose
          )
        } else {
          boot_obj <- bootMI(X_rand,
            X_nons,
            weights,
            y_nons,
            family_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights_rand,
            mu_hat,
            svydesign,
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
        # mu_hat <- boot_obj$mu
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(
          nonprob = NA,
          prob = NA
        ))))
      } else {
        stop("Invalid method for variance estimation.")
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

    X <- rbind(X_nons, X_rand) # joint model matrix
    # if (is.null(pop_size))
    pop_size <- N_est_rand # estimated pop_size

    output[[k]] <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = SE))))
    OutcomeList[[k]]$method <- method_outcome
    if (control_inference$vars_selection == TRUE) {
      OutcomeList[[k]]$cve <- cve_outcome
    } else {
      NULL
    }
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f
  names(OutcomeList) <- outcomes$f
  names(pop_size) <- "pop_size"
  names(ys) <- all.vars(outcome_init[[2]])

  boot_sample <- if (control_inference$var_method == "bootstrap" & control_inference$keep_boot) {
    stat
  } else {
    NULL
  }
  if (!is.null(boot_sample) & is.matrix(boot_sample)) colnames(boot_sample) <- names(ys)

  structure(
    list(
      X = if (isTRUE(x)) X else NULL,
      y = if (isTRUE(y)) ys else NULL,
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
      outcome = OutcomeList,
      boot_sample = boot_sample
    ),
    class = c("nonprobsvy", "nonprobsvy_mi")
  )
}
