#' @useDynLib nonprobsvy
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats qnorm
#' @importFrom stats as.formula
#' @importFrom stats terms
#' @import RcppArmadillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp

nonprobIPW <- function(selection,
                       target,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_selection,
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection = controlSel(),
                       control_inference = controlInf(),
                       start,
                       verbose,
                       x,
                       y,
                       se,
                       ...){

  h <- control_selection$h
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  var_method <- control_inference$var_method
  est_method <- control_selection$est_method_sel
  num_boot <- control_inference$num_boot
  var_selection <- control_inference$vars_selection
  lambda_control <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  nfolds <- control_selection$nfolds
  eps <- control_selection$epsilon
  if(!(target[3] == "NULL()")) stop("ill-defined formula for the target")
  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))
  outcomes <- ff(outcome)
  output <- list()
  ys <- list()
  if (se) {
    confidence_interval <- list()
    SE_values <- list()
  } else {
    confidence_interval <- NULL
    SE_values <- NULL
  }
  # formula for outcome variable if outcome defined
  # dependents <- paste(selection, collapse = " ")
  # outcome <- stats::as.formula(paste(outcome[2], dependents))

    if (is.null(pop_totals) && !is.null(svydesign)) {

      model <- model_frame(formula = outcomes$outcome[[1]],
                           data = data,
                           svydesign = svydesign)
      #y_nons <- model$y_nons

      X_nons <- model$X_nons
      X_rand <- model$X_rand
      nons_names <- model$nons_names
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

      if (var_selection == TRUE) {
        #X_stand <- cbind(1, ncvreg::std(X)) # standardization of variables before fitting
        X_stand <- ncvreg::std(X) # penalizing without an intercept
        prior_weights <- c(weights_rand, weights)

        method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
        method <- get_method(method_selection_function)
        inv_link <- method$make_link_inv

        n <- nrow(X)
        p <- ncol(X)
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
                                 lambda = lambda_control,
                                 pop_totals = pop_totals,
                                 verbose = verbose)
        min <- cv$min
        lambda <- cv$lambda
        theta_selected <- cv$theta_selected
        cve_selection <- cv$cv_error

        idx <- c(1, theta_selected + 2) # intercept plus selected variables
        #psel <- length(idx)
        X <- as.matrix(X[, idx])
        #X_design <- cbind(1, Xsel)
        X_nons <- X[loc_nons, ,drop = FALSE]
        X_rand <- X[loc_rand, ,drop = FALSE]
      }
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
                                        start = start,
                                        varcov = TRUE,
                                        control_selection = control_selection)

        estimation_method <- get_method(est_method)
        selection_model <- estimation_method$estimation_model(model = model_sel,
                                                        method_selection = method_selection)
        theta_hat <- selection_model$theta_hat
        #grad <- est_method_obj$grad
        hess <- selection_model$hess
        var_cov1 <- selection_model$var_cov1
        var_cov2 <- selection_model$var_cov2
        ps_nons <- selection_model$ps_nons
        est_ps_rand <- selection_model$est_ps_rand
        ps_nons_der <- selection_model$ps_nons_der
        est_ps_rand_der <- selection_model$est_ps_rand_der
        theta_standard_errors <- sqrt(diag(selection_model$variance_covariance))

        names(theta_hat) <- names(selection_model$theta_hat) <- colnames(X)
        weights_nons <- 1/ps_nons
        N <- sum(weights * weights_nons)

        # if (!is.null(pop_size)) {
        #   N <- pop_size
        # } else {
        #   N <- sum(weights * weights_nons)
        # }
      } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

        if (var_selection == FALSE) {
          if (!is.null(pop_totals)) pop_size <- pop_totals[1]
        }
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
        model <- model_frame(formula = outcomes$outcome[[1]],
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
        X <- rbind(model$X_rand, model$X_nons) # joint matrix

        if (var_selection == TRUE) {
          X_stand <- ncvreg::std(X) # penalizing without an intercept

          method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
          method <- get_method(method_selection_function)
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
                                   lambda = lambda_control,
                                   pop_totals = pop_totals[-1],
                                   verbose = verbose)
          min <- cv$min
          lambda <- cv$lambda
          theta_selected <- cv$theta_selected
          cve_selection <- cv$cv_error

          idx <-  c(1, theta_selected + 2) # intercept + selected variables
          psel <- length(idx)
          X <- as.matrix(X[, idx, drop = FALSE])
          X_nons <- X_nons[, idx, drop = FALSE]
          # data <- data[, idx, drop=FALSE]
          #colnames(X_design) <- c("(Intercept)", colnames(Xsel))

          ################ ESTIMATION
          pop_totals <- model$pop_totals[idx]
        }

        if (is.null(start)) {
          if (control_selection$start_type == "glm") {
            start <- start_fit(X = X, # <--- does not work with pop_totals
                               R = R,
                               weights = weights,
                               weights_rand = weights_rand,
                               method_selection = method_selection)
          } else if (control_selection$start_type == "naive") {
            start_h <- suppressWarnings(theta_h_estimation(R = R,
                                          X = X[,1,drop=FALSE],
                                          weights_rand = weights_rand,
                                          weights = weights,
                                          h = h,
                                          method_selection = method_selection,
                                          start = 0,
                                          maxit = maxit,
                                          pop_totals = pop_totals[1])$theta_h)
            start <- c(start_h, rep(0, ncol(X) - 1))
          }
        }

        h_object <- theta_h_estimation(R = R,
                                       X = X,
                                       weights_rand = weights_rand,
                                       weights = weights,
                                       h = h,
                                       method_selection = method_selection,
                                       start = start,
                                       maxit = maxit,
                                       pop_totals = pop_totals) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation

        theta_hat <- h_object$theta_h
        hess <- h_object$hess
        grad <- h_object$grad
        names(theta_hat) <- colnames(X)
        method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
        method <- get_method(method_selection_function)
        inv_link <- method$make_link_inv
        dinv_link <- method$make_link_inv_der
        eta_nons <- theta_hat %*% t(X_nons)
        ps_nons <- inv_link(eta_nons)
        ps_nons_der <- dinv_link(eta_nons)
        variance_covariance <- solve(-hess)
        theta_standard_errors <- sqrt(diag(variance_covariance))
        var_cov1 <- method$variance_covariance1
        var_cov2 <- method$variance_covariance2
        df_residual <- nrow(X_nons) - length(theta_hat)
        weights_nons <- 1/ps_nons
        N <- sum(weights * weights_nons)

        selection_model <- list(theta_hat = theta_hat,
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
        #######################################
      } else {
        stop("Please, provide svydesign object or pop_totals/pop_means.")
      }

  for (k in 1:outcomes$l) {
    if (is.null(pop_totals)) {
      y_nons <- model_frame(formula = outcomes$outcome[[k]],
                            data = data,
                            svydesign = svydesign,
                            pop_size = pop_size,
                            weights = weights,
                            flag = FALSE)$y_nons
    } else {
      y_nons <- model_frame(formula = outcomes$outcome[[k]],
                            data = data,
                            pop_totals = pop_totals,
                            pop_size = pop_size,
                            weights = weights,
                            flag = FALSE)$y_nons
    }
    ys$k <- as.numeric(y_nons) # TODO name to change
    mu_hat <- mu_hatIPW(y = y_nons,
                        weights = weights,
                        weights_nons = weights_nons,
                        N = ifelse(is.null(pop_size), N, pop_size)) # IPW estimator # consider using weighted.mean function
    #mu_hat <- weighted.mean(y_nons, w = weights * weights_nons)
    if (se) {
      if (var_method == "analytic") {
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
      } else if (var_method == "bootstrap") {

        if (control_inference$cores > 1) {
          boot_obj <- bootIPW_multicore(X_rand = X_rand,
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
                                        control_inference = control_inference,
                                        cores = control_inference$cores,
                                        verbose = verbose)
        } else {
          boot_obj <- bootIPW(X_rand = X_rand,
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
                             control_inference = control_inference,
                             verbose = verbose)
        }
        var <- boot_obj$var
        mu_hat <- boot_obj$mu
        SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      } else {
        stop("Invalid method for variance estimation.")
      }
      SE <- sqrt(var)
      alpha <- control_inference$alpha
      z <- stats::qnorm(1-alpha/2)
      # confidence interval based on the normal approximation
      confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * SE,
                                                                       upper_bound = mu_hat + z * SE
      ))))
    } else {
      SE <- NA
      confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = NA,
                                                                       upper_bound = NA
      ))))
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = NA, prob = NA))))
    }

    output[[k]] <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = SE))))
  }

  X <- rbind(X_nons, X_rand) # joint model matrix
  parameters <- matrix(c(theta_hat, theta_standard_errors),
                       ncol = 2,
                       dimnames = list(names(theta_hat),
                                       c("Estimate", "Std. Error")))
  weights_summary <- summary(as.vector(weights_nons))
  prop_scores <- c(ps_nons, est_ps_rand)
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f
  if (is.null(pop_size)) pop_size <- N # estimated pop_size
  names(pop_size) <- "pop_size"

  SelectionList <- list(coefficients = selection_model$theta_hat,
                        std_err = theta_standard_errors,
                        residuals = selection_model$residuals,
                        fitted_values = prop_scores,
                        link = selection_model$method,
                        linear_predictors = selection_model$eta,
                        aic = selection_model$aic,
                        weights = as.vector(weights_nons),
                        prior.weights = weights,
                        formula = selection,
                        df_residual = selection_model$df_residual,
                        log_likelihood = selection_model$log_likelihood)


  structure(
    list(X = if(isTRUE(x)) X else NULL,
         y = if(isTRUE(y)) ys else NULL,
         prob = prop_scores,
         weights = as.vector(weights_nons),
         control = list(control_selection = control_selection,
                        control_inference = control_inference),
         output = output,
         SE = SE_values,
         confidence_interval = confidence_interval,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         selection = SelectionList
  ),
  class = c("nonprobsvy", "nonprobsvy_ipw"))
}


mu_hatIPW <- function(y,
                      weights,
                      weights_nons,
                      N) {

  mu_hat <- 1/N * sum(weights * weights_nons * y)
  mu_hat

}


