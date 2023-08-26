#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats qnorm
#' @importFrom stats as.formula
#' @importFrom stats terms
#' @export
#' @rdname main_doc


nonprobIPW <- function(selection,
                       target,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size = NULL,
                       overlap,
                       method_selection,
                       family_selection = "binomial",
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection = controlSel(),
                       control_inference = controlInf(),
                       start,
                       verbose,
                       contrasts,
                       model,
                       x,
                       y,
                       ...){

  h <- control_selection$h
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  var_method <- control_inference$var_method
  est_method <- control_selection$est_method_sel
  num_boot <- control_inference$num_boot
  if(!(target[3] == "NULL()")) stop("ill-defined formula for the target")
  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))
  outcomes <- ff(outcome)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()

  # formula for outcome variable if outcome defined
  # dependents <- paste(selection, collapse = " ")
  # outcome <- stats::as.formula(paste(outcome[2], dependents))

  for (k in 1:outcomes$l) {
    if (is.null(pop_totals) && !is.null(svydesign)) {

      model <- model_frame(formula = outcomes$outcome[[k]],
                           data = data,
                           svydesign = svydesign)
      X_nons <- model$X_nons
      X_rand <- model$X_rand
      nons_names <- model$nons_names
      y_nons <- model$y_nons
      X <- rbind(X_rand, X_nons)

      R_nons <- rep(1, nrow(X_nons))
      R_rand <- rep(0, nrow(X_rand))
      R <- c(R_rand, R_nons)

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(X_nons)
      n_rand <- nrow(X_rand)

      #if (all(svydesign$prob) == 1) { # TODO
      #  if (!is.null(pop_size)) {
      #      ps_rand <- rep(sum(svydesign$prob)/pop_size, length(svydesign$prob))
      #    } else {
      #      ps_rand <- svydesign$prob/sum(svydesign$prob)
      #    }
      #} else {
      #  ps_rand <- svydesign$prob
      #}
      ps_rand <- svydesign$prob
      weights_rand <- 1/ps_rand
      weights_sum <- sum(weights_rand, weights)

      # case when samples overlap - to finish
      if (overlap) { # TODO
        overlap_idx_nons <- which(data[,control_selection$key] == 1)
        overlap_idx_rand <- which(svydesign$variables[,control_selection$key] == 1)
        overlap_model <- nonprobOv(X_nons,
                                   X_rand,
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
        N <- sum(weights_nons)
        mu_hat <-  mu_hatIPW(y = y_nons,
                             weights = weights,
                             weights_nons = weights_nons,
                             N = N)
        print(N)
        print(mu_hat)
        stop("model in development")
        # TODO objects
        selection <- NULL
        theta_hat <- overlap_model$parameters[,1]
        theta_standard_errors <- overlap_model$parameters[,2]
      } else {
        model_sel <- internal_selection(X = X,
                                        X_nons = X_nons,
                                        X_rand = X_rand,
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
        theta_hat <- selection$theta_hat
        #grad <- est_method_obj$grad
        hess <- selection$hess
        var_cov1 <- selection$var_cov1
        var_cov2 <- selection$var_cov2
        ps_nons <- selection$ps_nons
        est_ps_rand <- selection$est_ps_rand
        ps_nons_der <- selection$ps_nons_der
        est_ps_rand_der <- selection$est_ps_rand_der
        theta_standard_errors <- sqrt(diag(selection$variance_covariance))

        names(theta_hat) <- colnames(X)
        weights_nons <- 1/ps_nons

        if (!is.null(pop_size)) {
          N <- pop_size
        } else {
          N <- sum(weights * weights_nons)
        }

        mu_hat <- mu_hatIPW(y = y_nons,
                            weights = weights,
                            weights_nons = weights_nons,
                            N = N) # IPW estimator # consider using weighted.mean function
        #mu_hat <- weighted.mean(y_nons, w = weights * weights_rand)

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

        # names_pop <- names(pop_totals)
        # #pop_totals <- as.vector(pop_totals)
        # names(pop_totals) <- names_pop

        # model for outcome formula
        model <- model_frame(formula = outcomes$outcome[[k]],
                             data = data,
                             pop_totals = pop_totals)

        X_nons <- model$X_nons
        nons_names <- model$nons_names
        y_nons <- model$y_nons
        R <- rep(1, nrow(X_nons))
        n_nons <- nrow(X_nons)
        pop_totals <- model$pop_totals

        X_rand <- NULL
        est_ps_rand <- NULL
        est_ps_rand_der <- NULL
        ps_rand <- NULL
        n_rand <- 0
        weights_rand <- NULL
        log_likelihood <- "NULL"

        h_object <- theta_h_estimation(R = R,
                                       X = X_nons,
                                       weights_rand = weights_rand,
                                       weights = weights,
                                       h = h,
                                       method_selection = method_selection,
                                       maxit = maxit,
                                       pop_totals = pop_totals) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation

        theta_hat <- h_object$theta_h
        hess <- h_object$hess
        grad <- h_object$grad
        names(theta_hat) <- model$total_names
        method <- get_method(method_selection)
        inv_link <- method$make_link_inv
        dinv_link <- method$make_link_inv_der
        eta_nons <- theta_hat %*% t(X_nons)
        ps_nons <- inv_link(eta_nons)
        ps_nons_der <- dinv_link(eta_nons)
        N <- sum(weights * 1/ps_nons)
        variance_covariance <- solve(-hess)
        theta_standard_errors <- sqrt(diag(variance_covariance))
        var_cov1 <- method$variance_covariance1
        var_cov2 <- method$variance_covariance2
        df_residual <- nrow(X_nons) - length(theta_hat)
        weights_nons <- 1/ps_nons

        selection <- list(theta_hat = theta_hat,
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
                          log_likelihood = "NULL")

        mu_hat <- mu_hatIPW(model$y_nons, weights = weights, weights_nons = weights_nons, N = N)
      }
      else {
        stop("Please, provide svydesign object or pop_totals/pop_means.")
      }

      if (var_method == "analytic") {
        if (overlap) {
          var <- boot_overlap(X_rand = X_rand,
                              X_nons = X_nons,
                              y = y_nons,
                              weights_nons = weights_nons,
                              weights = weights,
                              mu_hat = mu_hat,
                              O_hat = O_hat,
                              L_hat = L_hat,
                              weights_rand = weights_rand,
                              weights_rnons = weights_rnons,
                              method_selection = method_selection,
                              dependency = control_selection$dependence,
                              N = N,
                              type = "IPW",
                              idx_nonprob = overlap_idx_nons,
                              idx_prob = overlap_idx_rand,
                              control = control_selection)
          var <- as.vector(var)
          SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
        } else {
          var_obj <- internal_varIPW(svydesign = svydesign,
                                     X_nons = X_nons,
                                     X_rand = X_rand,
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
                                     theta = theta_hat,
                                     h = h,
                                     var_cov1 = var_cov1,
                                     var_cov2 = var_cov2)

          var_nonprob <- var_obj$var_nonprob
          var_prob <- var_obj$var_prob
          var <- var_obj$var
          se_nonprob <- sqrt(var_nonprob)
          se_prob <- sqrt(var_prob)
          SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
        }
      } else if (var_method == "bootstrap") {

        if (control_inference$cores > 1) {
          var_obj <- bootIPW_multicore(X_rand = X_rand,
                                        X_nons = X_nons,
                                        y = y_nons,
                                        num_boot = num_boot,
                                        weights = weights,
                                        weights_rand = weights_rand,
                                        R = R,
                                        theta_hat = theta_hat,
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
                                        cores = control_inference$cores,
                                        verbose = verbose)
        } else {
          var_obj <- bootIPW(X_rand = X_rand,
                             X_nons = X_nons,
                             y = y_nons,
                             num_boot = num_boot,
                             weights = weights,
                             weights_rand = weights_rand,
                             R = R,
                             theta_hat = theta_hat,
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
        }

        var <- var_obj$boot_var
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      } else {
        stop("Invalid method for variance estimation.")
      }
    if (is.null(pop_size)) pop_size <- N # estimated pop_size
    se <- sqrt(var)

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                     upper_bound = mu_hat + z * se
                                                     ))))


    X <- rbind(X_nons, X_rand) # joint model matrix
    output[[k]] <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = se))))
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
                        control_inference = control_inference),
         output = output,
         SE = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         #log_likelihood = log_likelihood,
         #df_residual = df_residual,
         selection = selection
  ),
  class = c("nonprobsvy", "nonprobsvy_ipw"))
}


mu_hatIPW <- function(y,
                      weights,
                      weights_nons,
                      N) {

  mu_hat <- (1/N) * sum(weights * weights_nons * y)
  mu_hat

}

#' start_fit
#'
#' start_fit: Function for obtaining initial values for propensity score estimation
#'
#' @param X - a
#' @param R - a
#' @param weights - a
#' @param weights_rand - a
#' @param method_selection - a
#' @param control_selection - a

start_fit <- function(X,
                      R,
                      weights,
                      weights_rand,
                      method_selection,
                      control_selection = controlSel()) {

  weights_to_glm <- c(weights_rand, weights)

  start_model <- stats::glm.fit(x = X, #glm model for initial values in propensity score estimation
                                y = R,
                                weights = weights_to_glm, # to fix
                                #family = binomial(link = method_selection),
                                control = list(control_selection$epsilon,
                                               control_selection$maxit,
                                               control_selection$trace)

  )
  start_model$coefficients
}


