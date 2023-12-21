#' @useDynLib nonprobsvy
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats binomial
#' @importFrom stats terms
#' @importFrom ncvreg cv.ncvreg
#' @importFrom MASS ginv
#' @import Rcpp
#' @importFrom Rcpp evalCpp

nonprobDR <- function(selection,
                      outcome,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      method_selection,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_selection = controlSel(),
                      control_outcome = controlOut(),
                      control_inference = controlInf(),
                      start_outcome,
                      start_selection,
                      verbose,
                      x,
                      y,
                      se,
                      ...) {
  if (method_outcome %in% c("pmm", "nn")) {
    message("Bootstrap variance only, analytical version during implementation")
    control_inference$var_method <- "bootstrap"
  }
  outcome_init <- outcome
  h <- control_selection$h
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  est_method <- control_selection$est_method_sel
  var_method <- control_inference$var_method
  num_boot <- control_inference$num_boot
  bias_corr <- control_inference$bias_correction
  var_selection <- control_inference$vars_selection
  lambda_control <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  num_boot <- control_inference$num_boot
  eps <- control_selection$epsilon
  rep_type <- control_inference$rep_type
  ps_rand <- svydesign$prob
  weights_rand <- 1 / ps_rand

  outcomes <- ff(outcome)
  output <- list()
  OutcomeList <- list()
  ys <- list()
  if (se) {
    confidence_interval <- list()
    SE_values <- list()
  } else {
    confidence_interval <- NULL
    SE_values <- NULL
  }

  # Selection models
  if (is.null(pop_totals) && !is.null(svydesign)) {
    # TODO if (bias_corr == FALSE)
    # model for selection formula
    SelectionModel <- model_frame(
      formula = selection,
      data = data,
      svydesign = svydesign
    )
    # if (all(svydesign$prob) == 1) { # TODO
    #  if (!is.null(pop_size)) {
    #      ps_rand <- rep(sum(svydesign$prob)/pop_size, length(svydesign$prob))
    #    } else {
    #      ps_rand <- svydesign$prob/sum(svydesign$prob)
    #    }
    # } else {
    #  ps_rand <- svydesign$prob
    # }

    n_nons <- nrow(SelectionModel$X_nons)
    n_rand <- nrow(SelectionModel$X_rand)
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
    R_nons <- rep(1, n_nons)
    R_rand <- rep(0, n_rand)
    R <- c(R_rand, R_nons)
    loc_nons <- which(R == 1)
    loc_rand <- which(R == 0)
    weights_sum <- sum(weights_rand, weights)


    # y_rand <- vector(mode = "numeric", length = n_rand)
    # y <- c(y_rand, y_nons) # outcome variable for joint model
    n <- nrow(X)
    p <- ncol(X)
    N_rand <- sum(weights_rand)

    # Estimation for outcome model
    # model_out <- internal_outcome(outcome = outcome,
    #                               data = data,
    #                               weights = weights,
    #                               family_outcome = family_outcome)
    #
    # model_nons_coefs <- model_out$glm$coefficients
    # beta_statistics <- model_out$glm_summary$coefficients
    #
    # y_rand_pred <- as.numeric(OutcomeModel$X_rand %*% model_nons_coefs) # y_hat for probability sample
    # y_nons_pred <- model_out$glm$linear.predictors #as.numeric(X_nons %*% model_nons_coefs)

    # Estimation for selection model
    X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons) # joint model matrix
    ######  WORKING VERSION
    if (var_selection == TRUE) {
      X_stand <- ncvreg::std(X) # penalizing without an intercept
      prior_weights <- c(weights_rand, weights)

      method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
      method <- get_method(method_selection_function)
      inv_link <- method$make_link_inv

      # Cross-validation for variable selection
      if (verbose == TRUE & lambda_control == -1) {
        cat("Selection model\n")
      }
      cv <- cv_nonprobsvy_rcpp(
        X = X_stand,
        R = R,
        weights_X = prior_weights,
        method_selection = method_selection,
        h = h,
        maxit = maxit,
        eps = eps,
        lambda_min = lambda_min,
        nlambda = nlambda,
        nfolds = control_selection$nfolds,
        penalty = control_selection$penalty,
        a = switch(control_selection$penalty,
          SCAD = control_selection$a_SCAD,
          control_selection$a_MCP
        ),
        lambda = lambda_control,
        pop_totals = pop_totals,
        verbose = verbose
      )
      # theta_est <- cv$theta_est[cv$theta_est != 0]
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- c(0, cv$theta_selected + 1)
      cve_selection <- cv$cv_error
      lambda_selection <- cv$lambdas
      # names(theta_est) <- colnames(X)[theta_selected + 1]
      nlambda <- control_outcome$nlambda
    }
  } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {
    if (!is.null(pop_totals)) pop_size <- pop_totals[1]
    if (!is.null(pop_means)) { # TO consider
      if (!is.null(pop_size)) {
        pop_totals <- c(pop_size, pop_size * pop_means)
        names(pop_totals) <- c("(Intercept)", names(pop_means))
      } else {
        stop("pop_size must be defined when estimating with pop_means.")
      }
    }

    # model for selection formula
    SelectionModel <- model_frame(formula = selection, data = data, pop_totals = pop_totals)
    n_nons <- nrow(SelectionModel$X_nons)
    n_rand <- 0
    R <- rep(1, n_nons)
    X <- rbind(SelectionModel$X_nons, SelectionModel$X_rand) # joint model matrix

    ############# WORKING VERSION
    if (var_selection == TRUE) {
      X_stand <- ncvreg::std(X) # penalizing without an intercept

      method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
      method <- get_method(method_selection_function)
      inv_link <- method$make_link_inv

      n <- nrow(X)
      p <- ncol(X)

      # Cross-validation for variable selection
      if (verbose == TRUE & lambda_control == -1) {
        cat("Selection model\n")
      }
      cv <- cv_nonprobsvy_rcpp(
        X = X_stand, # TODO TO FIX
        R = R,
        weights_X = weights,
        method_selection = method_selection,
        h = h,
        maxit = maxit,
        eps = eps,
        lambda_min = lambda_min,
        nlambda = nlambda,
        nfolds = control_selection$nfolds,
        penalty = control_selection$penalty,
        a = switch(control_selection$penalty,
          SCAD = control_selection$a_SCAD,
          control_selection$a_MCP
        ),
        lambda = lambda_control,
        pop_totals = pop_totals[-1],
        verbose = verbose
      )
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- c(0, cv$theta_selected + 1)
      cve_selection <- cv$cv_error
    }
    ############
  } else {
    stop("Please, provide only one of svydesign object or pop_totals/pop_means.")
  }
  for (k in 1:outcomes$l) {
    outcome <- outcomes$outcome[[k]]

    if (is.null(pop_totals) && !is.null(svydesign)) {
      if (bias_corr == TRUE) {
        if (is.character(family_outcome)) {
          family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
          family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
          family_nonprobsvy <- family_nonprobsvy()
        }
        # Extract the terms from outcome and selection
        terms_out <- attr(terms(outcome, data = data), "term.labels")
        terms_sel <- attr(terms(selection, data = data), "term.labels")
        #
        # # Combine the terms
        combined_terms <- union(terms_out, terms_sel)
        combined_formula <- as.formula(paste(outcome[2], paste(combined_terms, collapse = " + "), sep = " ~ "))

        Model <- model_frame(
          formula = combined_formula,
          data = data,
          svydesign = svydesign
        )
        OutcomeModel <- SelectionModel <- Model
        n_nons <- nrow(Model$X_nons)
        n_rand <- nrow(Model$X_rand)
        R_nons <- rep(1, nrow(Model$X_nons))
        R_rand <- rep(0, nrow(Model$X_rand))
        R <- c(R_rand, R_nons)
        y_rand <- vector(mode = "numeric", length = n_rand)
        y <- c(y_rand, Model$y_nons) # outcome variable for joint model
        X <- rbind(Model$X_rand, Model$X_nons)

        ############ working version
        if (var_selection == TRUE) {
          nlambda <- control_outcome$nlambda
          if (verbose == TRUE) {
            cat("Outcome model\n")
          }
          beta <- ncvreg::cv.ncvreg(
            X = X[loc_nons, -1],
            y = Model$y_nons,
            penalty = control_outcome$penalty,
            family = family_outcome,
            nlambda = nlambda,
            trace = verbose,
            nfolds = control_outcome$nfolds,
            gamma = switch(control_outcome$penalty,
              SCAD = control_outcome$a_SCAD,
              control_outcome$a_MCP
            )
          )

          beta_est <- beta$fit$beta[, beta$min]
          beta_selected <- which(abs(beta_est) != 0) - 1
          beta_est <- beta_est[beta$fit$beta[, beta$min] != 0]
          cve_outcome <- beta$cve
          lambda_outcome <- beta$lambda
          lambda_min_outcome <- beta$lambda.min

          idx <- sort(unique(c(beta_selected[-1], theta_selected[-1]))) # excluding intercepts
          psel <- length(idx)
          Xsel <- as.matrix(X[, idx + 1, drop = FALSE])
          X <- cbind(1, Xsel)
          colnames(X) <- c("(Intercept)", colnames(Xsel))
          Model$X_rand <- X[loc_rand, ]
          Model$X_nons <- X[loc_nons, ]
        }
        estimation_model <- mm(
          X = X,
          y = y,
          weights = weights,
          weights_rand = weights_rand,
          R = R,
          n_nons = n_nons,
          n_rand = n_rand,
          method_selection = method_selection,
          family = family_nonprobsvy,
          start_selection = start_selection,
          start_outcome = start_outcome
        )

        selection_model <- estimation_model$selection
        OutcomeList[[k]] <- outcome_model <- estimation_model$outcome

        theta_hat <- selection_model$theta_hat
        grad <- selection_model$grad
        hess <- selection_model$hess
        ps_nons <- selection_model$ps_nons
        est_ps_rand <- selection_model$est_ps_rand
        ps_nons_der <- selection_model$ps_nons_der
        est_ps_rand_der <- selection_model$est_ps_rand_der
        theta_standard_errors <- sqrt(diag(selection_model$variance_covariance))
        names(theta_hat) <- names(selection_model$theta_hat) <- colnames(X)

        y_rand_pred <- outcome_model$y_rand_pred
        y_nons_pred <- outcome_model$y_nons_pred
        beta <- outcome_model$coefficients
        beta_errors <- sqrt(diag(outcome_model$variance_covariance))
        beta_statistics <- data.frame(beta = beta, beta_errors = beta_errors)
        sigma <- outcome_model$sigma_rand

        OutcomeList[[k]][c("sigma_rand", "y_rand_pred", "y_nons_pred")] <- NULL

        weights_nons <- 1 / ps_nons
        N_nons <- sum(weights * weights_nons)
        N_rand <- sum(weights_rand)
        y_nons <- Model$y_nons

        if (is.null(pop_size)) pop_size <- N_nons
        mu_hat <- mu_hatDR(
          y = Model$y_nons,
          y_nons = y_nons_pred,
          y_rand = y_rand_pred,
          weights = weights,
          weights_nons = weights_nons,
          weights_rand = weights_rand,
          N_nons = N_nons,
          N_rand = N_rand
        )

        # updating probability sample by adding y_hat variable
        svydesign <- stats::update(svydesign,
          .y_hat_MI = y_rand_pred
        )
      } else {
        # model for outcome formula
        OutcomeModel <- model_frame(
          formula = outcome,
          data = data,
          svydesign = svydesign
        )
        n_nons <- nrow(OutcomeModel$X_nons)
        n_rand <- nrow(OutcomeModel$X_rand)
        R_nons <- rep(1, n_nons)
        R_rand <- rep(0, n_rand)
        R <- c(R_rand, R_nons)
        loc_nons <- which(R == 1)
        loc_rand <- which(R == 0)
        weights_sum <- sum(weights_rand, weights)

        ############ working version
        if (var_selection == TRUE) {
          nlambda <- control_outcome$nlambda
          if (verbose == TRUE) {
            cat("Outcome model\n")
          }
          beta <- ncvreg::cv.ncvreg(
            X = X[loc_nons, -1],
            y = OutcomeModel$y_nons,
            penalty = control_outcome$penalty,
            family = family_outcome,
            nlambda = nlambda,
            trace = verbose,
            nfolds = control_outcome$nfolds,
            gamma = switch(control_outcome$penalty,
              SCAD = control_outcome$a_SCAD,
              control_outcome$a_MCP
            )
          )

          beta_est <- beta$fit$beta[, beta$min]
          beta_selected <- which(abs(beta_est) != 0) - 1
          beta_est <- beta_est[beta$fit$beta[, beta$min] != 0]
          cve_outcome <- beta$cve
          lambda_outcome <- beta$lambda
          lambda_min_outcome <- beta$lambda.min

          idx <- sort(unique(c(beta_selected[-1], theta_selected[-1]))) # excluding intercepts
          psel <- length(idx)
          Xsel <- as.matrix(X[, idx + 1, drop = FALSE])
          X <- cbind(1, Xsel)
          colnames(X) <- c("(Intercept)", colnames(Xsel))
          OutcomeModel$X_nons <- SelectionModel$X_nons <- X[loc_nons, ]
          OutcomeModel$X_rand <- SelectionModel$X_rand <- X[loc_rand, ]
        }

        model_sel <- internal_selection(
          X = X,
          X_nons = SelectionModel$X_nons,
          X_rand = SelectionModel$X_rand,
          weights = weights,
          weights_rand = weights_rand,
          R = R,
          method_selection = method_selection,
          optim_method = optim_method,
          h = h,
          est_method = est_method,
          maxit = maxit,
          control_selection = control_selection,
          start = start_selection
        )

        estimation_method <- get_method(est_method)
        selection_model <- estimation_method$estimation_model(
          model = model_sel,
          method_selection = method_selection
        )
        theta_hat <- selection_model$theta_hat
        # grad <- est_method_obj$grad
        hess <- selection_model$hess
        ps_nons <- selection_model$ps_nons
        est_ps_rand <- selection_model$est_ps_rand
        ps_nons_der <- selection_model$ps_nons_der
        est_ps_rand_der <- selection_model$est_ps_rand_der
        theta_standard_errors <- sqrt(diag(selection_model$variance_covariance))
        # log_likelihood <- est_method_obj$log_likelihood
        # df_residual <- est_method_obj$df_residual

        names(selection_model$theta_hat) <- names(theta_hat) <- colnames(X)
        weights_nons <- 1 / ps_nons
        N_nons <- sum(weights * weights_nons)
        N_rand <- sum(weights_rand)

        #############

        method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
        MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
        model_obj <- MethodOutcome(
          outcome = outcome,
          data = data,
          weights = weights,
          family_outcome = family_outcome,
          start_outcome = start_outcome,
          X_nons = OutcomeModel$X_nons,
          y_nons = OutcomeModel$y_nons,
          X_rand = OutcomeModel$X_rand,
          control = control_outcome,
          n_nons = n_nons,
          n_rand = n_rand,
          model_frame = OutcomeModel$model_frame,
          vars_selection = control_inference$vars_selection,
          pop_totals = pop_totals
        )

        y_rand_pred <- model_obj$y_rand_pred
        y_nons_pred <- model_obj$y_nons_pred
        outcome_model <- model_obj$model
        # beta_statistics <- model_obj$parameters
        sigma <- NULL
        OutcomeList[[k]] <- model_obj$model
        y_nons <- OutcomeModel$y_nons

        mu_hat <- mu_hatDR(
          y = OutcomeModel$y_nons,
          y_nons = y_nons_pred,
          y_rand = y_rand_pred,
          weights = weights,
          weights_nons = weights_nons,
          weights_rand = weights_rand,
          N_nons = N_nons,
          N_rand = N_rand
        )

        # updating probability sample by adding y_hat variable
        svydesign <- stats::update(svydesign,
          .y_hat_MI = y_rand_pred
        )
      }
    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {
      # model for outcome formula
      OutcomeModel <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)
      X_nons <- OutcomeModel$X_nons
      X_rand <- OutcomeModel$X_rand
      nons_names <- OutcomeModel$nons_names
      y_nons <- OutcomeModel$y_nons

      if (var_selection == TRUE) {
        nlambda <- control_outcome$nlambda
        if (verbose == TRUE) {
          cat("Outcome model\n")
        }
        beta <- ncvreg::cv.ncvreg(
          X = X[loc_nons, -1, drop = FALSE],
          y = OutcomeModel$y_nons,
          penalty = control_outcome$penalty,
          family = family_outcome,
          nlambda = nlambda,
          nfolds = control_outcome$nfolds,
          trace = verbose,
          gamma = switch(control_outcome$penalty,
            SCAD = control_outcome$a_SCAD,
            control_outcome$a_MCP
          )
        )

        beta_est <- beta$fit$beta[, beta$min]
        beta_selected <- which(abs(beta_est) != 0) - 1
        beta_est <- beta_est[beta$fit$beta[, beta$min] != 0]
        cve_outcome <- beta$cve
        lambda_outcome <- beta$lambda
        lambda_min_outcome <- beta$lambda.min

        # Estimating theta, beta parameters using selected variables
        # beta_selected <- beta_selected[-1] - 1

        idx <- sort(unique(c(beta_selected[-1], theta_selected[-1]))) # excluding intercepts
        psel <- length(idx)
        Xsel <- as.matrix(X[, idx + 1, drop = FALSE])
        X <- cbind(1, Xsel)
        colnames(X) <- c("(Intercept)", colnames(Xsel))
        X_nons <- OutcomeModel$X_nons <- SelectionModel$X_nons <- X[loc_nons, ]
        SelectionModel$pop_totals <- c(SelectionModel$pop_totals[1], SelectionModel$pop_totals[idx + 1])
      }

      if (is.null(start_selection)) {
        if (control_selection$start_type == "glm") {
          start_selection <- start_fit(
            X = SelectionModel$X_nons, # <--- does not work with pop_totals
            R = R,
            weights = weights,
            weights_rand = weights_rand,
            method_selection = method_selection
          )
        } else if (control_selection$start_type == "naive") {
          start_h <- suppressWarnings(theta_h_estimation(
            R = R,
            X = SelectionModel$X_nons[, 1, drop = FALSE],
            weights_rand = weights_rand,
            weights = weights,
            h = h,
            method_selection = method_selection,
            start = 0,
            maxit = maxit,
            pop_totals = SelectionModel$pop_totals[1]
          )$theta_h)
          start_selection <- c(start_h, rep(0, ncol(X) - 1))
        }
      }

      h_object <- theta_h_estimation(
        R = R,
        X = SelectionModel$X_nons,
        weights_rand = NULL,
        weights = weights,
        h = h,
        method_selection = method_selection,
        start = start_selection,
        maxit = maxit,
        pop_totals = SelectionModel$pop_totals
      )
      theta_hat <- h_object$theta_h
      hess <- h_object$hess
      grad <- h_object$grad
      names(theta_hat) <- SelectionModel$total_names
      method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
      method <- get_method(method_selection_function)
      inv_link <- method$make_link_inv
      dinv_link <- method$make_link_inv_der
      eta_nons <- theta_hat %*% t(SelectionModel$X_nons)
      ps_nons <- inv_link(eta_nons)
      ps_nons_der <- dinv_link(eta_nons)
      weights_nons <- 1 / ps_nons
      N_nons <- sum(weights * weights_nons)
      variance_covariance <- solve(-hess)
      theta_standard_errors <- sqrt(diag(variance_covariance))
      df_residual <- nrow(SelectionModel$X_nons) - length(theta_hat)
      # if(is.null(pop_size)) pop_size <- N_nons
      est_ps_rand <- NULL
      est_ps_rand_der <- NULL
      log_likelihood <- NA
      weights_rand <- NULL
      sigma <- NULL
      residuals <- as.vector(rep(1, n_nons) - ps_nons)

      selection_model <- list(
        theta_hat = theta_hat,
        hess = hess,
        grad = grad,
        ps_nons = ps_nons,
        est_ps_rand = est_ps_rand,
        ps_nons_der = ps_nons_der,
        est_ps_rand_der = est_ps_rand_der,
        variance_covariance = variance_covariance,
        # var_cov1 = var_cov1,
        # var_cov2 = var_cov2,
        df_residual = df_residual,
        log_likelihood = log_likelihood,
        eta = eta_nons,
        residuals = residuals,
        method = method
      )

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

      mu_hat <- 1 / N_nons * sum((1 / ps_nons) * (weights * (OutcomeModel$y_nons - y_nons_pred))) + y_rand_pred
    } else {
      stop("Please, provide only one of svydesign object or pop_totals/pop_means.")
    }
    ys[[k]] <- as.numeric(y_nons)

    if (se) {
      if (var_method == "analytic") { # TODO add estimator variance with model containg pop_totals to internal_varDR function
        var_obj <- internal_varDR(
          OutcomeModel = OutcomeModel, # consider add selection argument instead of separate arguments for selection objects
          SelectionModel = SelectionModel,
          y_nons_pred = y_nons_pred,
          weights = weights,
          weights_rand = weights_rand,
          method_selection = method_selection,
          control_selection = control_selection,
          ps_nons = ps_nons,
          theta = theta_hat,
          hess = hess,
          ps_nons_der = ps_nons_der,
          est_ps_rand = est_ps_rand,
          y_rand_pred = y_rand_pred,
          N_nons = N_nons,
          est_ps_rand_der = est_ps_rand_der,
          svydesign = svydesign,
          est_method = est_method,
          h = h,
          pop_totals = pop_totals,
          sigma = sigma,
          bias_correction = bias_corr
        )

        var_prob <- var_obj$var_prob
        var_nonprob <- var_obj$var_nonprob

        var <- var_prob + var_nonprob
        se_prob <- sqrt(var_prob)
        se_nonprob <- sqrt(var_nonprob)
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
      } else if (var_method == "bootstrap") {
        if (control_inference$cores > 1) {
          boot_obj <- bootDR_multicore(
            outcome = outcome,
            data = data,
            svydesign = svydesign,
            SelectionModel = SelectionModel,
            OutcomeModel = OutcomeModel,
            family_outcome = family_outcome,
            method_outcome = method_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights = weights,
            weights_rand = weights_rand,
            R = R,
            theta_hat,
            mu_hat = mu_hat,
            method_selection = method_selection,
            control_selection = control_selection,
            start_selection = start_selection,
            control_outcome = control_outcome,
            control_inference = control_inference,
            n_nons = n_nons,
            n_rand = n_rand,
            optim_method = optim_method,
            est_method = est_method,
            h = h,
            maxit = maxit,
            pop_totals = pop_totals,
            pop_size = pop_size,
            pop_means = pop_means,
            bias_correction = bias_corr,
            cores = control_inference$cores,
            verbose = verbose
          )
        } else {
          boot_obj <- bootDR(
            outcome = outcome,
            data = data,
            svydesign = svydesign,
            SelectionModel = SelectionModel,
            OutcomeModel = OutcomeModel,
            family_outcome = family_outcome,
            method_outcome = method_outcome,
            start_outcome = start_outcome,
            num_boot = num_boot,
            weights = weights,
            weights_rand = weights_rand,
            R = R,
            theta_hat,
            mu_hat = mu_hat,
            method_selection = method_selection,
            control_selection = control_selection,
            start_selection = start_selection,
            control_inference = control_inference,
            control_outcome = control_outcome,
            n_nons = n_nons,
            n_rand = n_rand,
            optim_method = optim_method,
            est_method = est_method,
            h = h,
            maxit = maxit,
            pop_totals = pop_totals,
            pop_size = pop_size,
            pop_means = pop_means,
            bias_correction = bias_corr,
            verbose = verbose
          )
        }
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = NA, prob = NA))))
        var <- boot_obj$var
        # mu_hat <- boot_obj$mu
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
    output[[k]] <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = SE))))
    # OutcomeList$std_err <- ifelse(method_outcome == "glm", beta_statistics[,2], NULL)
    # names(OutcomeList$std_err) <- names(OutcomeList$coefficients)
    # parameters <- matrix(c(theta_hat, theta_standard_errors),
    #                      ncol = 2,
    #                      dimnames = list(names(theta_hat),
    #                                      c("Estimate", "Std. Error")))
    OutcomeList[[k]]$method <- method_outcome
  }
  weights_summary <- summary(as.vector(weights_nons))
  prop_scores <- c(ps_nons, est_ps_rand)
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f
  names(OutcomeList) <- outcomes$f
  if (is.null(pop_size)) pop_size <- N_nons
  names(pop_size) <- "pop_size"
  names(ys) <- all.vars(outcome_init[[2]])

  SelectionList <- list(
    coefficients = selection_model$theta_hat,
    std_err = theta_standard_errors,
    residuals = selection_model$residuals,
    variance = as.vector(selection_model$variance),
    variance_covariance = selection_model$variance_covariance,
    fitted_values = prop_scores,
    family = selection_model$method,
    linear_predictors = selection_model$eta,
    aic = selection_model$aic,
    weights = as.vector(weights_nons),
    prior.weights = weights,
    formula = selection,
    df_residual = selection_model$df_residual,
    log_likelihood = selection_model$log_likelihood
  )
  # df.null = selection_model$df_null
  # converged)
  # TODO add imputed y to the outcome or separately

  structure(
    list(
      X = if (isTRUE(x)) X else NULL,
      y = if (isTRUE(y)) ys else NULL,
      prob = prop_scores,
      weights = as.vector(weights_nons),
      control = list(
        control_selection = control_selection,
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
      selection = SelectionList,
      boot_sample = if (control_inference$var_method == "bootstrap" & control_inference$keep_boot)
        boot_obj$stat else NULL
    ),
    class = c("nonprobsvy", "nonprobsvy_dr")
  )
}

mu_hatDR <- function(y,
                     y_nons,
                     y_rand,
                     weights,
                     weights_nons,
                     weights_rand,
                     N_nons,
                     N_rand) {
  mu_hat <- 1 / N_nons * sum(weights * weights_nons * (y - y_nons)) + 1 / N_rand * sum(weights_rand * y_rand)
  mu_hat
}
