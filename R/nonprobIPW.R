#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprobIPW} fits model for propensity score inference based on non-probability surveys using various methods.
#'
#' \loadmathjax
#'
#' @param selection `formula`, the selection (propensity) equation.
#' @param target `formula` with target variables.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals an optional `named vector` with population totals.
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
  est_method <- control_selection$est_method_sel

  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents)) # TODO reformulate

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
                                    varcov = TRUE)

    estimation_method <- get_method(est_method)
    est_method_obj <- estimation_method$estimation_model(model = model_sel,
                                                         method_selection = method_selection)
    theta_hat <- est_method_obj$theta_hat
    grad <- est_method_obj$grad
    hess <- est_method_obj$hess
    var_cov1 <- est_method_obj$var_cov1
    var_cov2 <- est_method_obj$var_cov2
    ps_nons <- est_method_obj$ps_nons
    est_ps_rand <- est_method_obj$est_ps_rand
    ps_nons_der <- est_method_obj$ps_nons_der
    est_ps_rand_der <- est_method_obj$est_ps_rand_der
    theta_standard_errors <- sqrt(diag(est_method_obj$variance_covariance))
    log_likelihood <- est_method_obj$log_likelihood
    df_residual <- est_method_obj$df_residual

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

  } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

    if (!is.null(pop_means)) {
      pop_totals <- pop_size * pop_means
    }

    names_pop <- names(pop_totals)
    pop_totals <- as.vector(pop_totals)
    names(pop_totals) <- names_pop

    # model for outcome formula
    model <- model_frame(formula = outcome,
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
    ps_rand <- NULL
    n_rand <- NULL
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
    theta_standard_errors <- sqrt(diag(solve(-hess)))
    var_cov1 <- method$variance_covariance1
    var_cov2 <- method$variance_covariance2
    df_residual <- nrow(X_nons) - length(theta_hat)
    weights_nons <- 1/ps_nons

    mu_hat <- mu_hatIPW(model$y_nons, weights = weights, weights_nons = weights_nons, N = N)
  }
  else {
    stop("Please, provide svydesign object or pop_totals/pop_means.")
  }

  if (var_method == "analytic") {
    var_obj <- internal_varIPW(svydesign = svydesign,
                               selection_formula = selection,
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
    SE_values <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
  } else if (var_method == "bootstrap") {
    var_obj <- bootIPW(X_rand = X_rand,
                       X_nons = X_nons,
                       y = y_nons,
                       family_outcome = family_outcome,
                       num_boot = 500,
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
                       pop_totals = pop_totals)
    var <- var_obj$boot_var
    SE_values <- "not computed for bootstrap variance"
  } else {
    stop("Invalid method for variance estimation.")
  }

  # case when samples overlap - to finish
  if (control_selection$overlap) { # TODO

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

  if (is.null(pop_size)) pop_size <- N # estimated pop_size

  mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)
  se <- sqrt(var)

  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  confidence_interval <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                   upper_bound = mu_hat + z * se
                                                   ))))


  X <- rbind(X_nons, X_rand) # joint model matrix
  output <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = se))))
  parameters <- data.frame("Estimate" = theta_hat,
                           "Std. Error" = theta_standard_errors,
                            row.names = names(theta_hat))
  weights_summary <- summary(as.vector(weights_nons))
  prop_scores <- c(ps_nons, est_ps_rand)


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
         log_likelihood = log_likelihood,
         df_residual = df_residual
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
                                family = binomial(link = method_selection),
                                control = list(control_selection$epsilon,
                                               control_selection$maxit,
                                               control_selection$trace)
  )
  start_model$coefficients

}


