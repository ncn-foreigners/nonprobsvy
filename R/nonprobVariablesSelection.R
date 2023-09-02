#' @useDynLib nonprobsvy
#' @importFrom MASS ginv
#' @importFrom ncvreg cv.ncvreg
#' @importFrom nleqslv nleqslv
#' @importFrom stats qnorm
#' @importFrom stats delete.response
#' @importFrom survey svymean
#' @importFrom Matrix Matrix
#' @importFrom stats terms
#' @import RcppArmadillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @export
#' @rdname main_doc


nonprobSel <- function(selection,
                       outcome,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_selection,
                       method_outcome,
                       family_selection = "binomial",
                       family_outcome = "gaussian",
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection = controlSel(),
                       control_outcome = controlOut(),
                       control_inference = controlInf(),
                       start,
                       verbose,
                       contrasts,
                       model,
                       x,
                       y,
                       ...) {

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  eps <- control_selection$epsilon
  maxit <- control_selection$maxit
  h <- control_selection$h
  lambda <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  num_boot <- control_inference$num_boot
  est_method <- control_selection$est_method_sel
  optim_method <- control_selection$optim_method
  bias_corr <- control_inference$bias_correction

  if (control_selection$a_SCAD <= 2 || control_outcome$a_SCAD <= 2) stop("a_SCAD must be greater than 2 for SCAD penalty")
  if (control_selection$a_MCP <= 1 || control_outcome$a_MCP <= 1) stop("a_MCP must be greater than 1 for MCP penalty")

  # XY_nons <- model.frame(outcome, data)
  # X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  # nons_names <- attr(terms(outcome, data = data), "term.labels")
  # if (all(nons_names %in% colnames(svydesign$variables))) {
  #   xx <- paste("~", paste(nons_names, collapse = "+"))
  #   outcome_rand <- as.formula(paste(outcome[2], xx))
  #   X_rand <- model.matrix(delete.response(terms(outcome_rand)), svydesign$variables[, nons_names]) # bug if formula is y~. #matrix of probability sample with intercept X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #
  # } else {
  #   stop("variable names in data and svydesign do not match")
  # } # TODO with dot_check for factor variables

  # Extract the terms from outcome and selection
  terms_out <- attr(terms(outcome, data = data), "term.labels")
  terms_sel <- attr(terms(selection, data = data), "term.labels")
  #
  # # Combine the terms
  combined_terms <- union(terms_out, terms_sel)
  combined_formula <- as.formula(paste(outcome[2], paste(combined_terms, collapse = " + "), sep = " ~ "))
  outcomes <- ff(combined_formula)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()

  for (k in 1:outcomes$l) {
    outcome <- outcomes$outcome[[k]]
    # OutcomeModel <- model_frame(formula = outcome,
    #                             data = data,
    #                             svydesign = svydesign)
    #
    # SelectionModel <- model_frame(formula = selection,
    #                               data = data,
    #                               svydesign = svydesign)
    ########################################################

    if (is.null(pop_totals) && !is.null(svydesign)) {

    Model <- model_frame(formula = outcome,
                         data = data,
                         svydesign = svydesign,
                         pop_size = pop_size,
                         weights = weights)
    OutcomeModel <- SelectionModel <- Model

    y_nons <- Model$y_nons
    ps_rand <- svydesign$prob
    weights_rand <- 1/ps_rand
    weights_sum <- sum(weights_rand, weights)

    R_nons <- rep(1, nrow(Model$X_nons))
    R_rand <- rep(0, nrow(Model$X_rand))
    R <- c(R_rand, R_nons)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise

    loc_nons <- which(R == 1)
    loc_rand <- which(R == 0)

    n_nons <- nrow(Model$X_nons)
    n_rand <- nrow(Model$X_rand)
    X <- rbind(Model$X_rand, Model$X_nons) # joint matrix
    #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
    X_stand <- ncvreg::std(X) # penalizing without an intercept
    prior_weights <- c(weights_rand, weights)

    method <- get_method(method_selection)
    inv_link <- method$make_link_inv

    y_rand <- vector(mode = "numeric", length = n_rand)
    y <- c(y_rand, y_nons) # outcome variable for joint model
    n <- nrow(X)
    p <- ncol(X)
    N_rand <- sum(weights_rand)

    # Cross-validation for variable selection
    if (verbose == TRUE & lambda == -1) {
      cat("Selection model\n")
    }
    cv <- cv_nonprobsvy_rcpp(X = X_stand,
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
                             a = switch(control_selection$penalty, SCAD = control_selection$a_SCAD, control_selection$a_MCP),
                             lambda = lambda,
                             pop_totals = pop_totals,
                             verbose = verbose)
    #theta_est <- cv$theta_est[cv$theta_est != 0]
    min <- cv$min
    lambda <- cv$lambda
    theta_selected <- c(0, cv$theta_selected + 1)
    #names(theta_est) <- colnames(X)[theta_selected + 1]

    nlambda <- control_outcome$nlambda
    if (verbose == TRUE) {
      cat("Outcome model\n")
    }
    beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1],
                              y = y_nons,
                              penalty = control_outcome$penalty,
                              family = family_outcome,
                              nlambda = nlambda,
                              trace = verbose,
                              nfolds = control_outcome$nfolds,
                              gamma = switch(control_outcome$penalty, SCAD = control_outcome$a_SCAD, control_outcome$a_MCP))

    beta_est <- beta$fit$beta[,beta$min]
    beta_selected <- which(abs(beta_est) != 0) - 1
    beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

    # Estimating theta, beta parameters using selected variables
    # beta_selected <- beta_selected[-1] - 1

    idx <- sort(unique(c(beta_selected[-1], theta_selected[-1]))) # excluding intercepts
    psel <- length(idx)
    Xsel <- as.matrix(X[, idx + 1])
    X_design <- cbind(1, Xsel)
    colnames(X_design) <- c("(Intercept)", colnames(Xsel))
    par0 <- rep(0, 2*(psel + 1))

    ################################## ESTIMATION on selected variables
    if (bias_corr == TRUE) {

      #mm <- get_method(control_inference$est_method)
      estimation_model <- mm(X = X_design,
                             y = y,
                             weights = weights,
                             weights_rand = weights_rand,
                             R = R,
                             n_nons = n_nons,
                             n_rand = n_rand,
                             method_selection = method_selection,
                             family = family_nonprobsvy)
      selection <- estimation_model$selection
      outcome <- estimation_model$outcome

      theta <- selection$theta_hat
      grad <- selection$grad
      hess <- selection$hess
      ps_nons <- selection$ps_nons
      est_ps_rand <- selection$est_ps_rand
      ps_nons_der <- selection$ps_nons_der
      est_ps_rand_der <- selection$est_ps_rand_der
      theta_errors <- sqrt(diag(selection$variance_covariance))

      y_rand_pred <- outcome$y_rand_pred
      y_nons_pred <- outcome$y_nons_pred
      beta <- outcome$beta_hat
      beta_errors <- sqrt(diag(outcome$variance_covariance))
      beta_statistics <- data.frame(beta = beta, beta_errors = beta_errors)
      sigma <- outcome$family$var

      weights_nons <- 1/ps_nons
      N_nons <- sum(weights * weights_nons)
      N_rand <- sum(weights_rand)

      OutcomeModel$X_nons <- cbind(1, OutcomeModel$X_nons[,idx+1])
      SelectionModel$X_nons <- cbind(1, SelectionModel$X_nons[,idx+1])

      OutcomeModel$X_rand <- cbind(1, OutcomeModel$X_rand[,idx+1])
      SelectionModel$X_rand <- cbind(1, SelectionModel$X_rand[,idx+1])

    } else {
      # OUTCOME MODEL
      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(outcome = outcome,
                                 data = data,
                                 weights = weights,
                                 family_outcome = family_outcome,
                                 X_nons = X_design[loc_nons,],
                                 y_nons = OutcomeModel$y_nons,
                                 X_rand = X_design[loc_rand,],
                                 control = control_outcome,
                                 n_nons = n_nons,
                                 n_rand = n_rand,
                                 model_frame = OutcomeModel$model_frame[,idx+1],
                                 vars_selection = control_inference$vars_selection,
                                 pop_totals = pop_totals)

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      outcome <- model_obj$model
      beta_statistics <- model_obj$parameters
      # beta_errors <- model_obj$parameters[,2]
      # beta <- model_obj$parameters[,1]
      sigma <- NULL

      # SELECTION MODEL
      model_sel <- internal_selection(X = X_design,
                                      X_nons = X_design[loc_nons,],
                                      X_rand = X_design[loc_rand,],
                                      weights = weights,
                                      weights_rand = weights_rand,
                                      R = R,
                                      method_selection = method_selection,
                                      optim_method = optim_method,
                                      h = h,
                                      est_method = est_method,
                                      maxit = maxit,
                                      control_selection = control_selection)

      estimation_method <- get_method(est_method)
      selection <- estimation_method$estimation_model(model = model_sel,
                                                      method_selection = method_selection)
      theta <- selection$theta_hat
      #grad <- est_method_obj$grad
      hess <- selection$hess
      ps_nons <- selection$ps_nons
      est_ps_rand <- selection$est_ps_rand
      ps_nons_der <- selection$ps_nons_der
      est_ps_rand_der <- selection$est_ps_rand_der
      theta_errors <- sqrt(diag(selection$variance_covariance))
      #log_likelihood <- est_method_obj$log_likelihood
      #df_residual <- est_method_obj$df_residual

      names(theta) <- colnames(X_design)
      weights_nons <- 1/ps_nons
      N_nons <- sum(weights * weights_nons)
      N_rand <- sum(weights_rand)

      OutcomeModel$X_nons <- cbind(1, OutcomeModel$X_nons[,idx+1])
      SelectionModel$X_nons <- cbind(1, SelectionModel$X_nons[,idx+1])

      OutcomeModel$X_rand <- cbind(1, OutcomeModel$X_rand[,idx+1])
      SelectionModel$X_rand <- cbind(1, SelectionModel$X_rand[,idx+1])
    }

    if(is.null(pop_size)) pop_size <- N_nons
    mu_hat <- mu_hatDR(y = y_nons,
                       y_nons = y_nons_pred,
                       y_rand = y_rand_pred,
                       weights = weights,
                       weights_nons = weights_nons,
                       weights_rand = weights_rand,
                       N_nons = N_nons,
                       N_rand = N_rand) #DR estimator

    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) { # TODO
      if (!is.null(pop_means)) {
        if (!is.null(pop_size)) {
          pop_totals <- pop_size * pop_means
        } else {
          stop("pop_size must be defined when estimating with pop_means.")
        }
      }
      Model <- model_frame(formula = outcome,
                           data = data,
                           pop_totals = pop_totals,
                           pop_size = pop_size,
                           weights = weights)
      OutcomeModel <- SelectionModel <- Model

      ######## VARS SELECTION ###########3
      y_nons <- Model$y_nons

      R = rep(1, nrow(Model$X_nons))

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(Model$X_nons)
      n_rand <- nrow(Model$X_rand)
      X <- rbind(Model$X_rand, Model$X_nons) # joint matrix
      #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
      X_stand <- ncvreg::std(X) # penalizing without an intercept

      method <- get_method(method_selection)
      inv_link <- method$make_link_inv

      n <- nrow(X)
      p <- ncol(X)

      # Cross-validation for variable selection
      if (verbose == TRUE & lambda == -1) {
        cat("Selection model\n")
      }
      cv <- cv_nonprobsvy_rcpp(X = X_stand, # TODO TO FIX
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
                               a = switch(control_selection$penalty, SCAD = control_selection$a_SCAD, control_selection$a_MCP),
                               lambda = lambda,
                               pop_totals = pop_totals[-1],
                               verbose = verbose)
      #theta_est <- cv$theta_est[cv$theta_est != 0]
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- c(0, cv$theta_selected + 1)
      #names(theta_est) <- colnames(X)[theta_selected + 1]

      nlambda <- control_outcome$nlambda
      if (verbose == TRUE) {
        cat("Outcome model\n")
      }
      beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1],
                                y = y_nons,
                                penalty = control_outcome$penalty,
                                family = family_outcome,
                                nlambda = nlambda,
                                nfolds = control_outcome$nfolds,
                                trace = verbose,
                                gamma = switch(control_outcome$penalty, SCAD = control_outcome$a_SCAD, control_outcome$a_MCP))

      beta_est <- beta$fit$beta[,beta$min]
      beta_selected <- which(abs(beta_est) != 0) - 1
      beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

      # Estimating theta, beta parameters using selected variables
      # beta_selected <- beta_selected[-1] - 1

      idx <- sort(unique(c(beta_selected[-1], theta_selected[-1]))) # excluding intercepts
      psel <- length(idx)
      Xsel <- as.matrix(X[, idx + 1])
      X_design <- cbind(1, Xsel)
      colnames(X_design) <- c("(Intercept)", colnames(Xsel))
      par0 <- rep(0, 2*(psel + 1))

      ################ ESTIMATION
      pop_totals_sel <- c(Model$pop_totals[1], Model$pop_totals[idx+1])

      h_object <- theta_h_estimation(R = R,
                                     X = X_design,
                                     weights_rand = NULL,
                                     weights = weights,
                                     h = h,
                                     method_selection = method_selection,
                                     maxit = maxit,
                                     pop_totals = pop_totals_sel)
      theta <- h_object$theta_h
      hess <- h_object$hess
      grad <- h_object$grad
      names(theta) <- colnames(X_design)
      method <- get_method(method_selection)
      inv_link <- method$make_link_inv
      dinv_link <- method$make_link_inv_der
      n_nons <- nrow(Model$X_nons)
      ps_nons <- inv_link(theta %*% t(X_design))
      ps_nons_der <- dinv_link(theta %*% t(X_design))
      weights_nons <- 1/ps_nons
      N_nons <- sum(weights * weights_nons)
      variance_covariance <- solve(-hess)
      theta_errors <- sqrt(diag(variance_covariance))
      df_residual <- nrow(Model$X_nons) - length(theta)
      if(is.null(pop_size)) pop_size <- N_nons
      n_rand <- 0
      est_ps_rand <- NULL
      est_ps_rand_der <- NULL
      log_likelihood <- "NULL"
      weights_rand <- NULL

      selection <- list(theta_hat = theta,
                        hess = hess,
                        grad = grad,
                        ps_nons = ps_nons,
                        est_ps_rand = est_ps_rand,
                        ps_nons_der =  ps_nons_der,
                        est_ps_rand_der = est_ps_rand_der,
                        variance_covariance = variance_covariance,
                        #var_cov1 = var_cov1,
                        #var_cov2 = var_cov2,
                        df_residual = df_residual,
                        log_likelihood = log_likelihood)


      #model_out <- internal_outcome(X_nons = OutcomeModel$X_nons,
      #                              X_rand =  OutcomeModel$pop_totals, # <--- pop_size is an intercept in the model
      #                              y = OutcomeModel$y_nons,
      #                              weights = weights,
      #                              family_outcome = family_outcome)

      glm <- stats::glm.fit(x = X_design, y = Model$y_nons, weights = weights)
      glm_summary <- summary.glm(glm)
      outcome <- list(glm = glm,
                      glm_summary = glm_summary)

      model_nons_coefs <- outcome$glm$coefficients
      beta_statistics <- outcome$glm_summary$coefficients
      # beta_errors <- beta_statistics[,2]
      # beta <- beta_statistics[,1]
      new_data <- data.frame(as.list(pop_totals_sel))
      names(new_data) <- names(pop_totals_sel)
      #y_rand_pred <- predict.glm(outcome$glm, newdata = new_data, type = "response") # TODO
      y_nons_pred <- outcome$glm$fitted.values
      y_rand_pred <- t(model_nons_coefs) %*% pop_totals_sel
      sigma <- NULL

      OutcomeModel$X_nons <- cbind(1, OutcomeModel$X_nons[,idx+1])
      SelectionModel$X_nons <- cbind(1, SelectionModel$X_nons[,idx+1])
      pop_totals <- pop_totals_sel

      mu_hat <- 1/N_nons * sum((1/ps_nons)*(weights * (Model$y_nons - y_nons_pred))) + 1/pop_size * y_rand_pred
    } else {
      stop("Please, provide only one of svydesign object or pop_totals/pop_means.")
    }

    if (control_inference$var_method == "analytic") {
      var_obj <- internal_varDR(OutcomeModel = OutcomeModel, # consider add selection argument instead of separate arguments for selection objects
                                SelectionModel = SelectionModel,
                                y_nons_pred = y_nons_pred,
                                weights = weights,
                                weights_rand = weights_rand,
                                method_selection = method_selection,
                                control_selection = control_selection,
                                ps_nons = ps_nons,
                                theta = theta,
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
                                bias_correction = bias_corr)

      var_prob <- var_obj$var_prob
      var_nonprob <- var_obj$var_nonprob

      var <- var_prob + var_nonprob
      se_prob <- sqrt(var_prob)
      se_nonprob <- sqrt(var_nonprob)
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
    } else if (control_inference$var_method == "bootstrap") {
      if (control_inference$cores > 1) {
        var <- bootDR_multicore(SelectionModel = SelectionModel,
                                OutcomeModel = OutcomeModel,
                                family_outcome = family_outcome,
                                num_boot = num_boot,
                                weights = weights,
                                weights_rand = weights_rand,
                                R = R,
                                theta_hat = theta,
                                mu_hat = mu_hat,
                                method_selection = method_selection,
                                control_selection = control_selection,
                                n_nons = n_nons,
                                n_rand = n_rand,
                                optim_method = optim_method,
                                est_method = est_method,
                                h = h,
                                maxit = maxit,
                                pop_totals = pop_totals,
                                pop_size = pop_size,
                                pop_means = pop_means,
                                cores = control_inference$cores,
                                verbose = verbose)
      } else {
        var <- bootDR(SelectionModel = SelectionModel,
                      OutcomeModel = OutcomeModel,
                      family_outcome = family_outcome,
                      num_boot = num_boot,
                      weights = weights,
                      weights_rand = weights_rand,
                      R = R,
                      theta_hat = theta,
                      mu_hat = mu_hat,
                      method_selection = method_selection,
                      control_selection = control_selection,
                      n_nons = n_nons,
                      n_rand = n_rand,
                      optim_method = optim_method,
                      est_method = est_method,
                      h = h,
                      maxit = maxit,
                      pop_totals = pop_totals,
                      pop_size = pop_size,
                      pop_means = pop_means,
                      verbose = verbose)
      }
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      var <- var$boot_var
    }
    ps <- c(est_ps_rand, ps_nons)
    se <- sqrt(var) # standard error
    output[[k]] <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                                upper_bound = mu_hat + z * se))))

    parameters <- matrix(c(theta, theta_errors),
                         ncol = 2,
                         dimnames = list(names(theta),
                                         c("Estimate", "Std. Error")))

    # beta <- matrix(c(beta, beta_errors),
    #                ncol = 2,
    #                dimnames = list(names(theta),
    #                                c("Estimate", "Std. Error")))

    probabilities_summary <- summary(as.vector(ps_nons))
    prop_scores <- as.vector(ps)
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f

  structure(
    list(X = X,
         prop_scores = prop_scores,
         weights = as.vector(weights_nons),
         control = list(control_selection = control_selection,
                        control_outcome = control_outcome,
                        control_inference = control_inference),
         output = output,
         confidence_interval = confidence_interval,
         SE_values = SE_values,
         parameters = parameters,
         beta = beta_statistics,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = N_nons,
         #df_residual = df_residual,
         outcome = outcome,
         selection = selection
         #log_likelihood = "NULL"
         ),
    class = c("nonprobsvy", "nonprobsvy_dr"))
}

#' @rdname main_doc
nonprobSelM <- function(outcome, # TODO add pop_totals
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_outcome,
                        family_outcome,
                        subset,
                        strata,
                        weights,
                        na_action,
                        control_outcome,
                        control_inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...) {

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  # TODO for pop_totals
  # XY_nons <- model.frame(outcome, data)
  # X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
  # nons_names <- attr(terms(outcome, data = data), "term.labels")
  # if (all(nons_names %in% colnames(svydesign$variables))) {
  #   X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept
  # } else {
  #   stop("variable names in data and svydesign do not match")
  # }


  outcomes <- ff(outcome)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()
  num_boot <- control_inference$num_boot
  for (k in 1:outcomes$l) {
    if (is.null(pop_totals) && !is.null(svydesign)) {
      outcome <- outcomes$outcome[[k]]
      Model <- model_frame(formula = outcome,
                           data = data,
                           svydesign = svydesign)
      # y_nons <- XY_nons[,1]
      y_nons <- Model$y_nons
      ps_rand <- svydesign$prob
      weights_rand <- 1/ps_rand
      N_est_rand <- sum(weights_rand)

      R_nons <- rep(1, nrow(Model$X_nons))
      R_rand <- rep(0, nrow(Model$X_rand))
      R <- c(R_rand, R_nons)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(Model$X_nons)
      n_rand <- nrow(Model$X_rand)
      y_rand <- vector(mode = "numeric", length = n_rand)
      y <- c(y_rand, y_nons) # outcome variable for joint model
      X <- rbind(Model$X_rand, Model$X_nons) # joint matrix
      #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
      prior_weights <- c(weights_rand, weights)

      nlambda <- control_outcome$nlambda
      beta <- ncvreg::cv.ncvreg(X = X[loc_nons, -1],
                                y = y_nons,
                                penalty = control_outcome$penalty,
                                family = family_outcome,
                                trace = verbose,
                                nlambda = nlambda)

      beta_est <- beta$fit$beta[,beta$min]
      beta_selected <- which(abs(beta_est) != 0) - 1
      beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

      X_design <- as.matrix(X[, beta_selected + 1])
      # colnames(X_design) <- c("(Intercept)", colnames(Xsel))
      X_rand <- X_design[loc_rand, ]
      X_nons <- X_design[loc_nons, ]
      # X <- rbind(X_rand, X_nons) # joint model matrix

      # multiroot <- nleqslv::nleqslv(x = par0, # TODO add user-specified parameters to control functions
      #                               fn = u_beta_mi,
      #                               method = "Broyden", # TODO consider the method
      #                               global = "qline",
      #                               xscalm = "fixed",
      #                               jacobian = TRUE,
      #                               R = R,
      #                               X = Xsel,
      #                               y = y,
      #                               weights = prior_weights,
      #                               family_nonprobsvy = family_nonprobsvy)
      # beta_sel <- beta_est
      # names(beta_sel) <- c("(Intercept)", nons_names[beta_selected])
      # N_est_rand <- sum(weights_rand)
      # #y_rand_pred <- as.numeric(X_rand %*% beta_sel) # y_hat for probability sample # consider predict function
      # #y_nons_pred <- as.numeric(X_nons %*% beta_sel)
      # df_residual <- nrow(X) - length(beta_sel)
      #
      # eta <- as.numeric(X %*% beta_sel)
      # y_hat <- family_nonprobsvy$mu(eta)
      # y_rand_pred <- y_hat[loc_rand]
      # y_nons_pred <- y_hat[loc_nons]
      # sigma <- family_nonprobsvy$variance(mu = y_hat, y = y[loc_rand])
      # residuals <- family_nonprobsvy$residuals(mu = y_rand_pred, y = y[loc_rand])
      #
      # # variance-covariance matrix for outcome model
      # vcov_outcome <- solve(t(X) %*% (sigma * X))
      # beta_errors <- sqrt(diag(vcov_outcome))
      #
      # outcome <- list(beta_hat = beta_sel,
      #                 grad = NULL,
      #                 hess = NULL,
      #                 #variance_covariance = vcov(model_out$glm),
      #                 df_residual = df_residual,
      #                 log_likelihood = "NULL")
      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(outcome = outcome,
                                 data = data,
                                 weights = weights,
                                 family_outcome = family_outcome,
                                 X_nons = X_nons,
                                 y_nons = y_nons,
                                 X_rand = X_rand,
                                 control = control_outcome,
                                 n_nons = n_nons,
                                 n_rand = n_rand,
                                 model_frame = Model$model_frame[,beta_selected+1],
                                 vars_selection = control_inference$vars_selection,
                                 pop_totals = pop_totals)

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      outcome <- model_obj$model
      # beta_statistics <- model_obj$parameters
      # beta_errors <- model_obj$parameters[,2]
      # beta <- model_obj$parameters[,1]
      sigma <- NULL
      svydesign <- stats::update(svydesign,
                                 y_hat_MI = y_rand_pred)

      mu_hat <- weighted.mean(x = y_rand_pred, w = weights_rand)

      if (is.null(pop_size)) pop_size <- N_est_rand # estimated pop_size
    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {
      if (!is.null(pop_means)) { # TO consider
        if (!is.null(pop_size)) {
          pop_totals <- c(pop_size, pop_size * pop_means)
          names(pop_totals) <- c("(Intercept)", names(pop_means))
        } else {
          stop("pop_size must be defined when estimating with pop_means.")
        }
      }
      Model <- model_frame(formula = outcome,
                           data = data,
                           pop_totals = pop_totals)
      # y_nons <- XY_nons[,1]
      y_nons <- Model$y_nons

      R <- rep(1, nrow(Model$X_nons)) # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise
      loc_nons <- which(R == 1)
      n_nons <- nrow(Model$X_nons)
      n_rand <- NULL
      weights_rand <- NULL
      X_rand <- Model$X_rand

      nlambda <- control_outcome$nlambda
      beta <- ncvreg::cv.ncvreg(X = Model$X_nons[,-1],
                                y = y_nons,
                                penalty = control_outcome$penalty,
                                family = family_outcome,
                                trace = verbose,
                                nlambda = nlambda)

      beta_est <- beta$fit$beta[,beta$min]
      beta_selected <- which(abs(beta_est) != 0) - 1
      beta_est <- beta_est[beta$fit$beta[,beta$min] != 0]

      X_nons <- Model$X_nons[, beta_selected + 1]
      X <- X_nons

      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      pop_totals <- pop_totals[beta_selected+1]
      model_obj <- MethodOutcome(outcome = outcome,
                                 data = data,
                                 weights = weights,
                                 family_outcome = family_outcome,
                                 X_nons = X_nons,
                                 y_nons = y_nons,
                                 X_rand = NULL,
                                 control = control_outcome,
                                 n_nons = n_nons,
                                 n_rand = NULL,
                                 model_frame = NULL,
                                 vars_selection = control_inference$vars_selection,
                                 pop_totals = pop_totals)

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      outcome <- model_obj$model
      # beta_statistics <- model_obj$parameters
      # beta_errors <- model_obj$parameters[,2]
      # beta <- model_obj$parameters[,1]
      sigma <- NULL
      N_est_rand <- pop_totals[1]
      mu_hat <- ifelse(method_outcome == "glm", as.vector(y_rand_pred/N_est_rand), y_rand_pred)
    }

    if (control_inference$var_method == "analytic"){
      var_obj <- internal_varMI(svydesign = svydesign,
                                X_nons = X_nons, # X_nons
                                X_rand = X_rand, # X_rand
                                y = y_nons,
                                y_pred = y_nons_pred,
                                weights_rand = weights_rand,
                                method = method_outcome,
                                n_rand = n_rand,
                                n_nons = n_nons,
                                N = N_est_rand,
                                family = family_outcome,
                                parameters = model_obj$parameters,
                                pop_totals = pop_totals)

      var_nonprob <- var_obj$var_nonprob
      var_prob <- var_obj$var_prob

      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
      # variance
      var <- var_nonprob + var_prob

    } else if (control_inference$var_method == "bootstrap") {
      # bootstrap variance
      var <- bootMI(X_nons = X_nons, # X_nons
                    X_rand = X_rand, # X_rand
                    weights = weights,
                    y = y_nons,
                    family_outcome = family_outcome,
                    num_boot = num_boot,
                    weights_rand = weights_rand,
                    mu_hat = mu_hat,
                    svydesign = svydesign,
                    rep_type = control_inference$rep_type,
                    method = method_outcome,
                    control = control_outcome,
                    pop_totals = pop_totals,
                    verbose = verbose)
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
    }

    se <- sqrt(var)
    output[[k]] <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # parameters <- matrix(c(beta, beta_errors),
    #                ncol = 2,
    #                dimnames = list(names(beta),
    #                                c("Estimate", "Std. Error")))

    # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                                upper_bound = mu_hat + z * se))))
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f


  structure(
    list(X = X,
         control = list(control_outcome = control_outcome,
                        control_inference = control_inference,
                        method_outcome = method_outcome),
         output = output,
         confidence_interval = confidence_interval,
         SE_values = SE_values,
         parameters = model_obj$parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         outcome = outcome
    ),
    class = c("nonprobsvy", "nonprobsvy_mi"))
}

#' @rdname main_doc
nonprobSelP <- function(selection,
                        target,
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_selection,
                        method_outcome,
                        family_selection,
                        family_outcome,
                        subset,
                        strata,
                        weights,
                        na_action,
                        control_selection,
                        control_outcome,
                        control_inference,
                        start,
                        verbose,
                        contrasts,
                        model,
                        x,
                        y,
                        ...) { #TODO add pop_totals

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  eps <- control_selection$epsilon
  maxit <- control_selection$maxit
  h <- control_selection$h
  lambda <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  nfolds <- control_selection$nfolds
  optim_method <- control_selection$optim_method
  est_method <- control_selection$est_method_sel
  var_method <- control_inference$var_method
  num_boot <- control_inference$num_boot


  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))

  outcomes <- ff(outcome)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()

  for (k in 1:outcomes$l) {

    Model <- model_frame(formula = outcomes$outcome[[k]],
                         data = data,
                         svydesign = svydesign,
                         pop_totals = pop_totals,
                         pop_size = pop_size,
                         weights = weights)

    # XY_nons <- model.frame(outcome, data)
    # X_nons <- model.matrix(XY_nons, data) #matrix of nonprobability sample with intercept
    # nons_names <- attr(terms(outcome, data = data), "term.labels")
    # nons_names <- colnames(X_nons)
    # if (all(nons_names %in% colnames(svydesign$variables))) {
    #   X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #X_rand <- model.matrix(delete.response(terms(outcome)) , svydesign$variables) bug if formula is y~. #matrix of probability sample with intercept
    # } else {
    #   stop("variable names in data and svydesign do not match")
    # }

    # TODO for pop_totals
    # if (names(XY_nons[1]) %in% colnames(svydesign$variables)) {
    #   #X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #X_rand <- model.matrix(delete.response(terms(outcome)) , svydesign$variables) bug if formula is y~. #matrix of probability sample with intercept
    #   XY_rand <- model.frame(outcome, svydesign$variables)
    #   X_rand <- model.matrix(XY_rand, svydesign$variables)
    # } else {
    #   stop("variable names in data and svydesign do not match")
    # }

    if (is.null(pop_totals)) {

      Model <- model_frame(formula = outcomes$outcome[[k]],
                           data = data,
                           svydesign = svydesign,
                           pop_size = pop_size,
                           weights = weights)

      y_nons <- Model$y_nons
      ps_rand <- svydesign$prob
      weights_rand <- 1/ps_rand

      R_nons <- rep(1, nrow(Model$X_nons))
      R_rand <- rep(0, nrow(Model$X_rand))
      R <- c(R_rand, R_nons)  # a vector of the binary indicator of belonging to the nonprobability sample; 1 if the unit belongs, 0 otherwise

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(Model$X_nons)
      n_rand <- nrow(Model$X_rand)
      X <- rbind(Model$X_rand, Model$X_nons) # joint matrix
      #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
      X_stand <- ncvreg::std(X) # penalizing without an intercept
      prior_weights <- c(weights_rand, weights)

      method <- get_method(method_selection)
      inv_link <- method$make_link_inv

      y_rand <- vector(mode = "numeric", length = n_rand)
      y <- c(y_rand, y_nons) # outcome variable for joint model
      n <- nrow(X)
      p <- ncol(X)
      N_rand <- sum(weights_rand)

      # Cross-validation for variable selection
      cv <- cv_nonprobsvy_rcpp(X = X_stand,
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
                               a = switch(control_selection$penalty, SCAD = control_selection$a_SCAD, control_selection$a_MCP),
                               lambda = lambda,
                               pop_totals = pop_totals,
                               verbose = verbose)
      #theta_est <- cv$theta_est[cv$theta_est != 0]
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- cv$theta_selected
      #names(theta_est) <- colnames(X)[theta_selected + 1]

      idx <- c(1, theta_selected + 2) # intercept plus selected variables
      #psel <- length(idx)
      X_design <- as.matrix(X[, idx])
      #X_design <- cbind(1, Xsel)

      model_sel <- internal_selection(X = X_design,
                                      X_nons = X_design[loc_nons,],
                                      X_rand = X_design[loc_rand,],
                                      weights = weights,
                                      weights_rand = weights_rand,
                                      R = R,
                                      method_selection = method_selection,
                                      optim_method = optim_method,
                                      h = h,
                                      est_method = est_method,
                                      maxit = maxit,
                                      varcov = TRUE,
                                      control_selection = control_selection)
      estimation_method <- get_method(est_method)
      selection <- estimation_method$estimation_model(model = model_sel,
                                                      method_selection = method_selection)
      theta <- selection$theta_hat
      hess <- selection$hess
      var_cov1 <- selection$var_cov1
      var_cov2 <- selection$var_cov2
      ps_nons <- selection$ps_nons
      est_ps_rand <- selection$est_ps_rand
      ps_nons_der <- selection$ps_nons_der
      est_ps_rand_der <- selection$est_ps_rand_der
      theta_errors <- sqrt(diag(selection$variance_covariance))
      weights_nons <- 1/ps_nons

      if (!is.null(pop_size)) {
        N <- pop_size
      } else {
        pop_size <- N <- sum(weights * weights_nons)
      }

      mu_hat <- mu_hatIPW(y = y_nons,
                          weights = weights,
                          weights_nons = weights_nons,
                          N = N)
    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) { # TODO
      if (!is.null(pop_means)) {
        if (!is.null(pop_size)) {
          pop_totals <- pop_size * pop_means
        } else {
          stop("pop_size must be defined when estimating with pop_means.")
        }
      }

      Model <- model_frame(formula = outcomes$outcome[[k]],
                           data = data,
                           pop_totals = pop_totals,
                           pop_size = pop_size,
                           weights = weights)

      ######## VARS SELECTION ###########3
      y_nons <- Model$y_nons

      R = rep(1, nrow(Model$X_nons))

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(Model$X_nons)
      n_rand <- nrow(Model$X_rand)
      X <- rbind(Model$X_rand, Model$X_nons) # joint matrix
      #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
      X_stand <- ncvreg::std(X) # penalizing without an intercept

      method <- get_method(method_selection)
      inv_link <- method$make_link_inv

      n <- nrow(X)
      p <- ncol(X)

      # Cross-validation for variable selection
      cv <- cv_nonprobsvy_rcpp(X = X_stand, # TODO TO FIX
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
                               a = switch(control_selection$penalty, SCAD = control_selection$a_SCAD, control_selection$a_MCP),
                               lambda = lambda,
                               pop_totals = pop_totals[-1],
                               verbose = verbose)
      #theta_est <- cv$theta_est[cv$theta_est != 0]
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- cv$theta_selected
      #names(theta_est) <- colnames(X)[theta_selected + 1]

      idx <-  c(1, theta_selected + 2) # intercept + selected variables
      psel <- length(idx)
      X_design <- as.matrix(X[, idx])
      #colnames(X_design) <- c("(Intercept)", colnames(Xsel))

      ################ ESTIMATION
      pop_totals <- Model$pop_totals[idx]

      h_object <- theta_h_estimation(R = R,
                                     X = X_design,
                                     weights_rand = NULL,
                                     weights = weights,
                                     h = h,
                                     method_selection = method_selection,
                                     maxit = maxit,
                                     pop_totals = pop_totals)
      theta <- h_object$theta_h
      hess <- h_object$hess
      grad <- h_object$grad
      names(theta) <- colnames(X_design)
      method <- get_method(method_selection)
      inv_link <- method$make_link_inv
      dinv_link <- method$make_link_inv_der
      n_nons <- nrow(Model$X_nons)
      ps_nons <- inv_link(theta %*% t(X_design))
      ps_nons_der <- dinv_link(theta %*% t(X_design))
      weights_nons <- 1/ps_nons
      N_nons <- sum(weights * weights_nons)
      variance_covariance <- solve(-hess)
      theta_errors <- sqrt(diag(variance_covariance))
      var_cov1 <- method$variance_covariance1
      var_cov2 <- method$variance_covariance2
      df_residual <- nrow(Model$X_nons) - length(theta)
      if(is.null(pop_size)) pop_size <- N_nons
      n_rand <- 0
      est_ps_rand <- NULL
      est_ps_rand_der <- NULL
      log_likelihood <- "NULL"
      weights_rand <- NULL

      selection <- list(theta_hat = theta,
                        hess = hess,
                        grad = grad,
                        ps_nons = ps_nons,
                        est_ps_rand = est_ps_rand,
                        ps_nons_der =  ps_nons_der,
                        est_ps_rand_der = est_ps_rand_der,
                        variance_covariance = variance_covariance,
                        var_cov1 = var_cov1,
                        var_cov2 = var_cov2,
                        df_residual = df_residual,
                        log_likelihood = log_likelihood)

      mu_hat <- mu_hatIPW(y = y_nons,
                          weights = weights,
                          weights_nons = weights_nons,
                          N = pop_size)
    } else {
      stop("Please, provide only one of svydesign object or pop_totals/pop_means.")
    }
    if (var_method == "analytic") {
      var_obj <- internal_varIPW(svydesign = svydesign,
                                 X_nons = X_design[loc_nons,],
                                 X_rand = X_design[loc_rand,],
                                 y_nons = y_nons,
                                 weights = weights,
                                 ps_nons = ps_nons,
                                 mu_hat = mu_hat,
                                 hess = hess,
                                 ps_nons_der = ps_nons_der,
                                 N = N,
                                 est_ps_rand = est_ps_rand,
                                 ps_rand = ps_rand,
                                 est_ps_rand_der = est_ps_rand_der,
                                 n_rand = n_rand,
                                 pop_size = pop_size,
                                 pop_totals = pop_totals,
                                 method_selection = method_selection,
                                 est_method = est_method,
                                 theta = theta,
                                 h = h,
                                 var_cov1 = var_cov1,
                                 var_cov2 = var_cov2)

      var_nonprob <- var_obj$var_nonprob
      var_prob <- var_obj$var_prob
      var <- var_obj$var
      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
    } else if (var_method == "bootstrap") {
      var_obj <- bootIPW(X_rand = X_design[loc_rand,],
                         X_nons = X_design[loc_nons,],
                         y = y_nons,
                         num_boot = num_boot,
                         weights = weights,
                         weights_rand = weights_rand,
                         R = R,
                         theta_hat = theta,
                         mu_hat = mu_hat,
                         method_selection = method_selection,
                         n_nons = n_nons,
                         n_rand = n_rand,
                         optim_method = optim_method,
                         est_method = est_method,
                         h = h,
                         maxit = maxit,
                         pop_size = pop_size,
                         pop_totals = pop_totals,
                         control_selection = control_selection,
                         verbose = verbose)
      var <- var_obj$boot_var
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
    } else {
      stop("Invalid method for variance estimation.")
    }

    se <- sqrt(var)
    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)
      # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                                upper_bound = mu_hat + z * se
    ))))
    parameters <- matrix(c(theta, theta_errors),
                         ncol = 2,
                         dimnames = list(names(theta),
                                         c("Estimate", "Std. Error")))

    output[[k]] <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))
    prop_scores <- c(ps_nons, est_ps_rand)
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f

  structure(
    list(X = X_design,
         prop_scores = prop_scores,
         weights = as.vector(weights_nons),
         control = list(control_selection = control_selection,
                        control_inference = control_inference),
         output = output,
         SE = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         selection = selection
    ),
    class = c("nonprobsvy", "nonprobsvy_ipw"))
}
