#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprobIPW} fits model for propensity score inference based on non-probability surveys using various methods.
#' \loadmathjax
#' @param selection `formula`, the selection (propensity) equation.
#' @param target `formula` with target variables.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_total an optional `named vector` with population totals.
#' @param pop_means an optional `named vector` with population means.
#' @param pop_size an optional `double` with population size.
#' @param method_selection a `character` with method for propensity scores estimation
#' @param family_selection a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights.
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores.
#' @param control_inference a list indicating parameters to use in inference based on probability and non-probability samples, contains parameters such as estimation method or variance method.
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation.
#' @param verbose verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... Additional, optional arguments.
#'
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats qnorm
#' @importFrom stats as.formula
#' @importFrom stats terms
#' @export



nonprobIPW <- function(selection,
                       target,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size = NULL,
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


  h <- control_selection$h_x
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  var_method <- control_inference$var_method
  smooth <- control_selection$smooth
  #weights <- rep.int(1, nrow(data)) # to remove

  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))

  # formula for outcome variable if outcome defined
  # dependents <- paste(selection, collapse = " ")
  # outcome <- stats::as.formula(paste(outcome[2], dependents))

  if (is.null(pop_totals) && !is.null(svydesign)) {
    model <- model_frame(formula = outcome,
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
    X <- rbind(X_rand, X_nons)

    ps_rand <- svydesign$prob
    weights_rand <- 1/ps_rand

    model_sel <- internal_selection(X = X,
                                    X_nons = X_nons,
                                    X_rand = X_rand,
                                    weights = weights,
                                    weights_rand = weights_rand,
                                    R = R,
                                    method_selection = method_selection,
                                    optim_method = optim_method,
                                    h = h,
                                    smooth = smooth,
                                    maxit = maxit,
                                    varcov = TRUE)

    if (!smooth) {

      maxLik_nons_obj <- model_sel$maxLik_nons_obj
      maxLik_rand_obj <- model_sel$maxLik_rand_obj
      log_likelihood <- model_sel$log_likelihood # maximum of the loglikelihood function
      theta_hat <- model_sel$theta

      ps_nons <- maxLik_nons_obj$ps
      est_ps_rand <- maxLik_rand_obj$ps
      hess <- maxLik_nons_obj$hess
      var_cov1 <- model_sel$var_cov1
      var_cov2 <- model_sel$var_cov2


      if (method_selection == "probit") { # for probit model, propensity score derivative is required
        ps_nons_der <- maxLik_nons_obj$psd
        est_ps_rand_der <- maxLik_rand_obj$psd
      }

    } else {
      theta_hat <- model_sel$theta_hat
      hess <- model_sel$hess
      grad <- model_sel$grad
      ps_nons <- model_sel$ps_nons
      est_ps_rand <- model_sel$est_ps_rand
      ps_nons_der <- model_sel$ps_nons_der
      est_ps_rand_der <- model_sel$est_ps_rand_der
      var_method <- "bootstrap"
      #TO DO - variance estimation for theta_h
    }

    names(theta_hat) <- c("(Intercept)", nons_names)
    weights_nons <- 1/ps_nons

    if (!is.null(pop_size)) {
      N <- pop_size
    } else {
      N <- sum(weights_nons)
    }

    mu_hat <- mu_hatIPW(y = y_nons,
                        weights = weights_nons,
                        N = N) # IPW estimator

    if (var_method == "analytic") {
     var_obj <- internal_varIPW(X_nons = X_nons,
                                X_rand = X_rand,
                                y_nons = y_nons,
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
                                method_selection = method_selection,
                                var_cov1 = var_cov1,
                                var_cov2 = var_cov2)

      var_nonprob <- var_obj$var_nonprob
      var_prob <- var_obj$var_prob
      var <- var_obj$var
      theta_hat_var <- var_obj$theta_hat_var
      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)
    } else if (var_method == "bootstrap") {
      var <- bootIPW(X_rand = X_rand,
                     X_nons = X_nons,
                     y = y_nons,
                     family_outcome = family_outcome,
                     num_boot = 500,
                     weights = weights,
                     weights_rand = weights_rand,
                     R = R,
                     mu_hat = mu_hat,
                     method_selection = method_selection,
                     n_nons = n_nons,
                     n_rand = n_rand,
                     optim_method = optim_method,
                     smooth = smooth,
                     h = h,
                     maxit = maxit,
                     pop_size = pop_size,
      )
      inf <- "not computed for bootstrap variance"
    } else {
      stop("Invalid method for variance estimation.")
    }

  } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

    if (!is.null(pop_totals)) {
      pop_totals <- pop_size * pop_means
    }

    # model for outcome formula
    model <- model_frame(formula = outcome,
                         data = data,
                         pop_totals = pop_totals)

    h_object <- theta_h_estimation(R = R,
                                   X = X_sel,
                                   weights_rand = weights_rand,
                                   weights = weights,
                                   h = h,
                                   method_selection = method_selection,
                                   maxit = maxit) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation

    theta_hat <- h_object$theta_h
    hess_h <- h_object$hess
    grad_h <- h_object$grad
    names(theta_hat) <- c("(Intercept)", model$nons_names)
    method <- get_method(method_selection)
    inv_link <- method$make_link_inv
    ps_nons <- inv_link(theta_hat %*% t(model$X_nons))
    N_nons <- sum(1/ps_nons)

    mu_hat <- mu_hatIPW(model$y_nons, weights = 1/ps_nons, N = N_nons)
    var <- 0
    se_nonprob <- 0
    se_prob <- 0
  }
  else {
    stop("Please, provide svydesign object or pop_totals/pop_means.")
  }

  # case when samples overlap - to finish
  if (control_selection$overlap) {

    weights <- nonprobOv(X_nons,
                         X_rand,
                         weights_rand,
                         dependent = control_selection$dependence,
                         method_selection)$weights

    O_hat <- nonprobOv(X_nons,
                       X_rand,
                       weights_rand,
                       dependent = control_selection$dependence,
                       method_selection)$O_hat

    L_hat <- nonprobOv(X_nons,
                       X_rand,
                       weights_rand,
                       dependent = control_selection$dependence,
                       method_selection)$L_hat

    weights_rnons <- nonprobOv(X_nons,
                         X_rand,
                         weights_rand,
                         dependent = control_selection$dependence,
                         method_selection)$weights_rnons

    N <- sum(weights)

    mu_hat_Ov <-  mu_hatIPW(y = y_nons,
                            weights = weights,
                            N = N)

    var <- boot_overlap(X_rand = X_rand,
                        X_nons = X_nons,
                        y = y_nons,
                        weights = weights,
                        O_hat = O_hat,
                        L_hat = L_hat,
                        weights_rand = weights_rand,
                        weights_rnons = weights_rnons,
                        dependency = control_selection$dependence,
                        N = N)
    var <- as.vector(var)
  }

  mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)
  se <- sqrt(var)

  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  ci <- c(mu_hat - z * se, mu_hat + z * se)

  structure(
    list(mean = mu_hat,
         #VAR = var,
         SE = se,
         #VAR_nonprob = V1,
         #VAR_prob = V2,
         SE_nonprob = ifelse(var_method == "analytic", se_nonprob, inf),
         SE_prob = ifelse(var_method == "analytic", se_prob, inf),
         #variance_covariance = V_mx,
         CI = ci,
         theta = theta_hat
         #theta_variance = theta_hat_var,
         #pearson_residuals = pearson_residuals,
         #deviance_residuals = deviance_residuals,
         #log_likelihood = log_likelihood
  ),
  class = "Inverse probability weighted")

}


mu_hatIPW <- function(y,
                      weights,
                      N) {

  mu_hat <- (1/N) * sum(weights * y)
  mu_hat

}

#' start_fit
#'
#' start_fit: Function for obtaining initial values for propensity score estimation
#'
#' @param X - a
#' @param R - a
#' @param weights - a
#' @param d - a
#' @param method_selection - a
#' @param control_selection - a

start_fit <- function(X,
                      R,
                      weights,
                      d,
                      method_selection,
                      control_selection = controlSel()) {

  weights <- c(weights, d)

  start_model <- stats::glm.fit(x = X, #glm model for initial values in propensity score estimation
                                y = R,
                                #weights = c(weights, d), # to fix
                                family = binomial(link = method_selection),
                                control = list(control_selection$epsilon,
                                               control_selection$maxit,
                                               control_selection$trace)
  )
  start <- start_model$coefficients
  start
}


