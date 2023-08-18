#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats binomial
#' @importFrom stats terms
#' @importFrom MASS ginv
#' @export
#' @rdname main_doc


nonprobDR <- function(selection,
                      outcome,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      overlap,
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

  h <- control_selection$h_x
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  est_method <- control_selection$est_method_sel
  var_method <- control_inference$var_method
  num_boot <- control_inference$num_boot
  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand

  outcomes <- ff(outcome)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()

  for (k in 1:outcomes$l) {
    outcome <- outcomes$outcome[[k]]
    if (is.null(pop_totals) && !is.null(svydesign)) {

      if (overlap) { # TODO
        overlap_idx_nons <- which(data[,control_selection$key] == 1)
        overlap_idx_rand <- which(svydesign$variables[,control_selection$key] == 1)
        # model for outcome formula
        OutcomeModel <- model_frame(formula = outcome,
                                    data = data,
                                    svydesign = svydesign)
        #model for selection formula
        SelectionModel <- model_frame(formula = selection,
                                      data = data,
                                      svydesign = svydesign)
        X <- rbind(SelectionModel$X_nons, SelectionModel$X_rand) # joint model matrix
        N_rand <- sum(weights_rand)

        n_nons <- nrow(OutcomeModel$X_nons)
        n_rand <- nrow(OutcomeModel$X_rand)
        R_nons <- rep(1, n_nons)
        R_rand <- rep(0, n_rand)
        R <- c(R_rand, R_nons)
        loc_nons <- which(R == 1)
        loc_rand <- which(R == 0)
        weights_sum <- sum(weights_rand, weights)
        overlap_model <- nonprobOv(SelectionModel$X_nons,
                                   SelectionModel$X_rand,
                                   weights = weights,
                                   weights_rand,
                                   dependent = control_selection$dependence,
                                   method_selection,
                                   key_var_prob = svydesign$variables[,control_selection$key],
                                   idx_nonprob = overlap_idx_nons,
                                   idx_prob = overlap_idx_rand,
                                   control = control_selection)

        O_hat <- overlap_model$O_hat
        ps_nons <- overlap_model$ps_nons
        est_ps_rand <- overlap_model$est_ps_rand
        weights_nons <- 1/ps_nons
        L_hat <- overlap_model$L_hat
        weights_rnons <- overlap_model$weights_rnons
        N_nons <- sum(weights * weights_nons)
        MethodOutcome <- get(method_outcome, mode = "function", envir = parent.frame())
        model_obj <- MethodOutcome(outcome = outcome,
                                   data = data,
                                   weights = weights,
                                   family_outcome = family_outcome,
                                   X_nons = OutcomeModel$X_nons,
                                   y_nons = OutcomeModel$y_nons,
                                   X_rand = OutcomeModel$X_rand,
                                   control = control_outcome,
                                   n_nons = n_nons,
                                   n_rand = n_rand,
                                   model_frame = OutcomeModel$model_frame,
                                   vars_selection = control_inference$vars_selection)
        # TODO objects
        Selection <- NULL
        outcome <- NULL
        theta_hat <- overlap_model$parameters[,1]
        theta_standard_errors <- overlap_model$parameters[,2]

        y_rand_pred <- model_obj$y_rand_pred
        y_nons_pred <- model_obj$y_nons_pred
        outcome <- model_obj$model
        beta_statistics <- model_obj$parameters
        if(is.null(pop_size)) pop_size <- N_nons

        mu_hat <- mu_hatDR(y = OutcomeModel$y_nons,
                           y_nons = y_nons_pred,
                           y_rand = y_rand_pred,
                           weights = weights,
                           weights_nons = weights_nons,
                           weights_rand = weights_rand,
                           N_nons = N_nons,
                           N_rand = N_rand)
        print(mu_hat)
        stop("model in development")
      } else {

        if (est_method == "mm") {
          # TODO
          if(is.character(family_outcome)) {
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

          Model <- model_frame(formula = combined_formula,
                               data = data,
                               svydesign = svydesign)
          OutcomeModel <- SelectionModel <- Model
          n_nons <- nrow(Model$X_nons)
          n_rand <- nrow(Model$X_rand)
          R_nons <- rep(1, nrow(Model$X_nons))
          R_rand <- rep(0, nrow(Model$X_rand))
          R <- c(R_rand, R_nons)
          y_rand <- vector(mode = "numeric", length = n_rand)
          y <- c(y_rand, Model$y_nons) # outcome variable for joint model
          X <- rbind(Model$X_rand, Model$X_nons)


          mm <- get_method(control_inference$est_method)
          estimation_model <- mm(X = X,
                                 y = y,
                                 weights = weights,
                                 weights_rand = weights_rand,
                                 R = R,
                                 n_nons = n_nons,
                                 n_rand = n_rand,
                                 method_selection = method_selection,
                                 family = family_nonprobsvy)

          Selection <- estimation_model$selection
          outcome <- estimation_model$outcome

          theta_hat <- Selection$theta_hat
          grad <- Selection$grad
          hess <- Selection$hess
          ps_nons <- Selection$ps_nons
          est_ps_rand <- Selection$est_ps_rand
          ps_nons_der <- Selection$ps_nons_der
          est_ps_rand_der <- Selection$est_ps_rand_der
          theta_standard_errors <- sqrt(diag(Selection$variance_covariance))

          y_rand_pred <- outcome$y_rand_pred
          y_nons_pred <- outcome$y_nons_pred
          beta <- outcome$beta_hat
          beta_errors <- sqrt(diag(outcome$variance_covariance))
          beta_statistics <- data.frame(beta = beta, beta_errors = beta_errors)
          sigma <- outcome$family$var

          weights_nons <- 1/ps_nons
          N_nons <- sum(weights * weights_nons)
          N_rand <- sum(weights_rand)

          if(is.null(pop_size)) pop_size <- N_nons
          mu_hat <- mu_hatDR(y = Model$y_nons,
                             y_nons = y_nons_pred,
                             y_rand = y_rand_pred,
                             weights = weights,
                             weights_nons = weights_nons,
                             weights_rand = weights_rand,
                             N_nons = N_nons,
                             N_rand = N_rand) # DR estimator - consider using weighted.mean function

          # updating probability sample by adding y_hat variable
          svydesign <- stats::update(svydesign,
                                     .y_hat_MI = y_rand_pred)
        } else {
          # model for outcome formula
          OutcomeModel <- model_frame(formula = outcome,
                                      data = data,
                                      svydesign = svydesign)
          #model for selection formula
          SelectionModel <- model_frame(formula = selection,
                                        data = data,
                                        svydesign = svydesign)
          #if (all(svydesign$prob) == 1) { # TODO
          #  if (!is.null(pop_size)) {
          #      ps_rand <- rep(sum(svydesign$prob)/pop_size, length(svydesign$prob))
          #    } else {
          #      ps_rand <- svydesign$prob/sum(svydesign$prob)
          #    }
          #} else {
          #  ps_rand <- svydesign$prob
          #}

          n_nons <- nrow(OutcomeModel$X_nons)
          n_rand <- nrow(OutcomeModel$X_rand)
          R_nons <- rep(1, n_nons)
          R_rand <- rep(0, n_rand)
          R <- c(R_rand, R_nons)
          loc_nons <- which(R == 1)
          loc_rand <- which(R == 0)
          weights_sum <- sum(weights_rand, weights)

          MethodOutcome <- get(method_outcome, mode = "function", envir = parent.frame())
          model_obj <- MethodOutcome(outcome = outcome,
                                     data = data,
                                     weights = weights,
                                     family_outcome = family_outcome,
                                     X_nons = OutcomeModel$X_nons,
                                     y_nons = OutcomeModel$y_nons,
                                     X_rand = OutcomeModel$X_rand,
                                     control = control_outcome,
                                     n_nons = n_nons,
                                     n_rand = n_rand,
                                     model_frame = OutcomeModel$model_frame,
                                     vars_selection = control_inference$vars_selection,
                                     pop_totals = pop_totals)

          y_rand_pred <- model_obj$y_rand_pred
          y_nons_pred <- model_obj$y_nons_pred
          outcome <- model_obj$model
          beta_statistics <- model_obj$parameters
          sigma <- NULL

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
          X_nons <- SelectionModel$X_nons
          X_rand <- SelectionModel$X_rand
          X <- rbind(SelectionModel$X_rand, SelectionModel$X_nons) # joint model matrix

          model_sel <- internal_selection(X = X,
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
                                          control_selection = control_selection)

          estimation_method <- get_method(est_method)
          Selection <- estimation_method$estimation_model(model = model_sel,
                                                          method_selection = method_selection)
          theta_hat <- Selection$theta_hat
          #grad <- est_method_obj$grad
          hess <- Selection$hess
          ps_nons <- Selection$ps_nons
          est_ps_rand <- Selection$est_ps_rand
          ps_nons_der <- Selection$ps_nons_der
          est_ps_rand_der <- Selection$est_ps_rand_der
          theta_standard_errors <- sqrt(diag(Selection$variance_covariance))
          #log_likelihood <- est_method_obj$log_likelihood
          #df_residual <- est_method_obj$df_residual

          names(theta_hat) <- colnames(X)
          weights_nons <- 1/ps_nons
          N_nons <- sum(weights * weights_nons)
          N_rand <- sum(weights_rand)

          if(is.null(pop_size)) pop_size <- N_nons
          mu_hat <- mu_hatDR(y = OutcomeModel$y_nons,
                             y_nons = y_nons_pred,
                             y_rand = y_rand_pred,
                             weights = weights,
                             weights_nons = weights_nons,
                             weights_rand = weights_rand,
                             N_nons = N_nons,
                             N_rand = N_rand) # DR estimator - consider using weighted.mean function

          # updating probability sample by adding y_hat variable
          svydesign <- stats::update(svydesign,
                                     .y_hat_MI = y_rand_pred)
        }
      }
      } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

        if (!is.null(pop_means)) { # TO consider
          if (!is.null(pop_size)) {
            pop_totals <- c(pop_size, pop_size * pop_means)
            names(pop_totals) <- c("(Intercept)", names(pop_means))
          } else {
            stop("pop_size must be defined when estimating with pop_means.")
          }
        }

        # model for outcome formula
        OutcomeModel <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)
        #model for selection formula
        SelectionModel <- model_frame(formula = selection, data = data, pop_totals = pop_totals)
        X <- rbind(SelectionModel$X_nons, SelectionModel$X_rand) # joint model matrix

        h_object <- theta_h_estimation(R = rep(1, nrow(SelectionModel$X_nons)),
                                       X = SelectionModel$X_nons,
                                       weights_rand = NULL,
                                       weights = weights,
                                       h = h,
                                       method_selection = method_selection,
                                       maxit = maxit,
                                       pop_totals = SelectionModel$pop_totals)
        theta_hat <- h_object$theta_h
        hess <- h_object$hess
        grad <- h_object$grad
        names(theta_hat) <- SelectionModel$total_names
        method <- get_method(method_selection)
        inv_link <- method$make_link_inv
        dinv_link <- method$make_link_inv_der
        n_nons <- nrow(SelectionModel$X_nons)
        ps_nons <- inv_link(theta_hat %*% t(SelectionModel$X_nons))
        ps_nons_der <- dinv_link(theta_hat %*% t(SelectionModel$X_nons))
        weights_nons <- 1/ps_nons
        N_nons <- sum(weights * weights_nons)
        variance_covariance <- solve(-hess)
        theta_standard_errors <- sqrt(diag(variance_covariance))
        df_residual <- nrow(SelectionModel$X_nons) - length(theta_hat)
        if(is.null(pop_size)) pop_size <- N_nons
        n_rand <- 0
        est_ps_rand <- NULL
        est_ps_rand_der <- NULL
        log_likelihood <- "NULL"
        weights_rand <- NULL
        sigma <- NULL

        Selection <- list(theta_hat = theta_hat,
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

        if(is.character(family_outcome)) {
          family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
          family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
          family_nonprobsvy <- family_nonprobsvy()
        }

        outcome <- internal_outcome(outcome = outcome,
                                    data = data,
                                    weights = weights,
                                    family_outcome = family_outcome)

        model_nons_coefs <- outcome$glm$coefficients
        beta_statistics <- outcome$glm_summary$coefficients

        eta <- OutcomeModel$pop_totals %*% model_nons_coefs
        y_rand_pred <- family_nonprobsvy$mu(eta)
        y_nons_pred <- outcome$glm$fitted.values

        # new_data <- data.frame(as.list(OutcomeModel$pop_totals))
        # names(new_data) <- names(OutcomeModel$pop_totals)
        # y_rand_pred <- predict.glm(outcome$glm, newdata = new_data, type = "response")
        # y_rand_pred <- t(model_nons_coefs) %*% OutcomeModel$pop_totals

        mu_hat <- 1/N_nons * sum((1/ps_nons)*(weights * (OutcomeModel$y_nons - y_nons_pred))) + 1/pop_size * y_rand_pred
      } else {
        stop("Please, provide only one of svydesign object or pop_totals/pop_means.")
      }

    if (var_method == "analytic") { # TODO add estimator variance with model containg pop_totals to internal_varDR function
      if (overlap) {
        var <- boot_overlap(X_rand = SelectionModel$X_rand,
                            X_nons = SelectionModel$X_nons,
                            y = OutcomeModel$y_nons,
                            weights_nons = weights_nons,
                            weights = weights,
                            mu_hat = mu_hat,
                            O_hat = O_hat,
                            L_hat = L_hat,
                            weights_rand = weights_rand,
                            weights_rnons = weights_rnons,
                            method_selection = method_selection,
                            family_outcome = family_outcome,
                            dependency = control_selection$dependence,
                            N = N_nons,
                            type = "DR",
                            idx_nonprob = overlap_idx_nons,
                            idx_prob = overlap_idx_rand,
                            control = control_selection)
        var <- as.vector(var)
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      } else {
        var_obj <- internal_varDR(OutcomeModel = OutcomeModel, # consider add selection argument instead of separate arguments for selection objects
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
                                  sigma = sigma) # TODO

        var_prob <- var_obj$var_prob
        var_nonprob <- var_obj$var_nonprob

        var <- var_prob + var_nonprob
        se_prob <- sqrt(var_prob)
        se_nonprob <- sqrt(var_nonprob)
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
      }
    } else if (var_method == "bootstrap") {
      if (control_inference$cores > 1) {
        var_obj <- bootDR_multicore(SelectionModel = SelectionModel,
                                    OutcomeModel = OutcomeModel,
                                    family_outcome = family_outcome,
                                    num_boot = num_boot,
                                    weights = weights,
                                    weights_rand = weights_rand,
                                    R = R,
                                    theta_hat,
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
                                    cores = control_inference$cores)
      } else {
        var_obj <- bootDR(SelectionModel = SelectionModel,
                          OutcomeModel = OutcomeModel,
                          family_outcome = family_outcome,
                          num_boot = num_boot,
                          weights = weights,
                          weights_rand = weights_rand,
                          R = R,
                          theta_hat,
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
                          pop_means = pop_means)
      }
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      var <- var_obj$boot_var
    } else {
      stop("Invalid method for variance estimation.")
    }
    se <- sqrt(var)
    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)
    # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                                upper_bound = mu_hat + z * se))))
    output[[k]] <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))
    parameters <- matrix(c(theta_hat, theta_standard_errors),
                         ncol = 2,
                         dimnames = list(names(theta_hat),
                                         c("Estimate", "Std. Error")))
    weights_summary <- summary(as.vector(weights_nons))
    prop_scores <- c(ps_nons, est_ps_rand)
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
                        control_inference = control_inference,
                        method_outcome = method_outcome),
         output = output,
         SE_values = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         beta = beta_statistics,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         #log_likelihood = log_likelihood,
         #df_residual = df_residual,
         outcome = outcome, # TO change
         selection = Selection
         ),
    class = c("nonprobsvy", "nonprobsvy_dr"))
}

#' mu_hatDR
#
#' mu_hatDR: Function for outcome variable estimation based on doubly robust estimation
#'
#' @param y - a
#' @param y_nons - a
#' @param y_rand - a
#' @param weights - a
#' @param weights_nons - a
#' @param weights_rand - a
#' @param N_nons - a
#' @param N_rand - a

mu_hatDR <- function(y,
                     y_nons,
                     y_rand,
                     weights,
                     weights_nons,
                     weights_rand,
                     N_nons,
                     N_rand) {

  mu_hat <- 1/N_nons * sum(weights * weights_nons * (y - y_nons)) + 1/N_rand * sum(weights_rand * y_rand)
  mu_hat
}

