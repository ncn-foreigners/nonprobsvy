#' @import mathjaxr
NULL
#' @title Inference with the non-probability survey samples.
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#'
#' @description \code{nonprobDR} fits model for doubly-robust inference based on non-probability surveys using various methods.
#'
#' \loadmathjax
#'
#' @param selection `formula`, the selection (propensity) equation.
#' @param outcome `formula`, the outcome equation.
#' @param data an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals an optional `named vector` with population totals.
#' @param pop_means an optional `named vector` with population means.
#' @param pop_size an optional `double` with population size.
#' @param method_selection a `character` with method for propensity scores estimation
#' @param method_outcome a `character` with method for response variable estimation
#' @param family_selection a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata an optional `vector` specifying strata.
#' @param weights an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... Additional, optional argumnents.
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats binomial
#' @importFrom stats terms
#' @importFrom MASS ginv
#' @export


nonprobDR <- function(selection,
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

  h <- control_selection$h_x
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  est_method <- control_selection$est_method_sel
  var_method <- control_inference$var_method


  if (is.null(pop_totals) && !is.null(svydesign)) {
    # model for outcome formula
    OutcomeModel <- model_frame(formula = outcome,
                                data = data,
                                svydesign = svydesign)


    #model for selection formula
    SelectionModel <- model_frame(formula = selection,
                                  data = data,
                                  svydesign = svydesign)
    X_sel <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
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

    n_nons <- nrow(OutcomeModel$X_nons)
    n_rand <- nrow(OutcomeModel$X_rand)
    R_nons <- rep(1, n_nons)
    R_rand <- rep(0, n_rand)
    R <- c(R_rand, R_nons)
    loc_nons <- which(R == 1)
    loc_rand <- which(R == 0)
    weights_rand <- 1/ps_rand

    # Estimation for outcome model
    model_out <- internal_outcome(outcome = outcome,
                                  data = data,
                                  weights = weights,
                                  family_outcome = family_outcome)

    model_nons_coefs <- model_out$glm$coefficients
    beta_statistics <- model_out$glm_summary$coefficients

    y_rand_pred <- as.numeric(OutcomeModel$X_rand %*% model_nons_coefs) # y_hat for probability sample
    y_nons_pred <- model_out$glm$linear.predictors #as.numeric(X_nons %*% model_nons_coefs)

    # Estimation for selection model
    X_nons <- SelectionModel$X_nons
    X_rand <- SelectionModel$X_rand
    model_sel <- internal_selection(X = X_sel,
                                    X_nons = SelectionModel$X_nons,
                                    X_rand = SelectionModel$X_rand,
                                    weights = weights,
                                    weights_rand = weights_rand,
                                    R = R,
                                    method_selection = method_selection,
                                    optim_method = optim_method,
                                    h = h,
                                    est_method = est_method,
                                    maxit = maxit)

    estimation_method <- get_method(est_method)
    est_method_obj <- estimation_method$estimation_model(model = model_sel,
                                                         method_selection = method_selection)
    theta_hat <- est_method_obj$theta_hat
    grad <- est_method_obj$grad
    hess <- est_method_obj$hess
    ps_nons <- est_method_obj$ps_nons
    est_ps_rand <- est_method_obj$est_ps_rand
    ps_nons_der <- est_method_obj$ps_nons_der
    est_ps_rand_der <- est_method_obj$est_ps_rand_der
    theta_standard_errors <- sqrt(diag(est_method_obj$variance_covariance))
    log_likelihood <- est_method_obj$log_likelihood
    df_residual <- est_method_obj$df_residual

    names(theta_hat) <- colnames(X_sel)
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
                       N_rand = N_rand) #DR estimator # consider using weighted.mean function

    # updating probability sample by adding y_hat variable
    svydesign <- stats::update(svydesign,
                               .y_hat_MI = y_rand_pred)

    if (var_method == "analytic") {

      var_obj <- internal_varDR(OutcomeModel = OutcomeModel,
                                SelectionModel = SelectionModel,
                                y_nons_pred = y_nons_pred,
                                weights = weights,
                                method_selection = method_selection,
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
                                h = h)

      var_prob <- var_obj$var_prob
      var_nonprob <- var_obj$var_nonprob

      var <- var_prob + var_nonprob
      se_prob <- sqrt(var_prob)
      se_nonprob <- sqrt(var_nonprob)
      SE_values <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))

    } else if (var_method == "bootstrap") {
      var_obj <- bootDR(SelectionModel = SelectionModel,
                        OutcomeModel = OutcomeModel,
                        family_outcome = family_outcome,
                        num_boot = 500,
                        weights = weights,
                        weights_rand = weights_rand,
                        R = R,
                        theta_hat,
                        mu_hat = mu_hat,
                        method_selection = method_selection,
                        n_nons = n_nons,
                        n_rand = n_rand,
                        optim_method = optim_method,
                        est_method = est_method,
                        h = h,
                        maxit = maxit
                        )
      SE_values <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
      var <- var_obj$boot_var
    } else {
      stop("Invalid method for variance estimation.")
    }

  } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

    # Consider variance for probability component
    if (!is.null(pop_means)) {
      pop_totals <- pop_size * pop_means
    }

    # model for outcome formula
    OutcomeModel <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)
    #model for selection formula
    SelectionModel <- model_frame(formula = selection, data = data, pop_totals = pop_totals)

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
    N_est <- sum(weights * weights_nons)
    variance_covariance <- solve(-hess)
    theta_standard_errors <- sqrt(diag(variance_covariance))
    df_residual <- nrow(SelectionModel$X_nons) - length(theta_hat)
    if(is.null(pop_size)) pop_size <- N_est
    n_rand <- NULL
    est_ps_rand <- NULL
    est_ps_rand_der <- NULL
    log_likelihood <- "NULL"

    model_sel <- list(theta_hat = theta_hat,
                      hess = hess,
                      grad = grad,
                      ps_nons = ps_nons,
                      est_ps_rand = est_ps_rand,
                      ps_nons_der =  ps_nons_der,
                      est_ps_rand_der = est_ps_rand_der,
                      variance_covariance = variance_covariance,
                      #var_cov1 = var_cov1,
                      #var_cov2 = var_cov2,
                      df_residual = df_residual)


    #model_out <- internal_outcome(X_nons = OutcomeModel$X_nons,
    #                              X_rand =  OutcomeModel$pop_totals, # <--- pop_size is an intercept in the model
    #                              y = OutcomeModel$y_nons,
    #                              weights = weights,
    #                              family_outcome = family_outcome)

    model_out <- internal_outcome(outcome = outcome,
                                  data = data,
                                  weights = weights,
                                  family_outcome = family_outcome)

    model_nons_coefs <- model_out$glm$coefficients
    beta_statistics <- model_out$glm_summary$coefficients
    y_rand_pred <- sum(OutcomeModel$pop_totals * model_nons_coefs) # consider %*% instead of sum()
    y_nons_pred <- model_out$glm$linear.predictors


    mu_hat <- 1/N_est * sum((1/ps_nons)*(weights * (OutcomeModel$y_nons - y_nons_pred))) + 1/pop_size * y_rand_pred

    # var_prob <- as.vector(attr(svydesign_mean, "var")) # to consider

    var_nonprob <- switch(method_selection,
                  "logit" = 1/N_est^2 * sum((1 - ps_nons)/ps_nons^2 * (weights*(OutcomeModel$y_nons - y_nons_pred))^2),
                  "cloglog" = 0,
                  "probit" = 0) # variance based on section 4.2 in the article

    #svydesign <- svydesign(ids = ~1, probs = 1, data = data.frame(y = y_rand_pred))
    #svydesign_mean <- survey::svymean(~y, svydesign)
    #var_prob <- as.vector(attr(svydesign_mean, "var")) # probability component
    var_prob <- 1/pop_size^2 * y_rand_pred
    se_prob <- sqrt(var_prob)
    se_nonprob <- sqrt(var_nonprob)
    var <- var_nonprob + var_prob
    SE_values <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))

  } else {
    stop("Please, provide svydesign object or pop_totals/pop_means.")
  }

  # case when samples overlap - TODO
  if (control_selection$overlap) {

    weights_nons <- nonprobOv(SelectionModel$X_nons,
                              SelectionModel$X_rand,
                              weights_rand,
                              dependent = control_selection$dependence,
                              method_selection)
    N <- sum(weights)

    mu_hat_Ov <-  mu_hatDR(y = OutcomeModel$y_nons,
                           y_nons = y_nons_pred,
                           y_rand = y_rand_pred,
                           weights_nons = weights,
                           weights_rand = weights_rand,
                           N_nons = N,
                           N_rand = N_est_rand)

  }

  X <- rbind(SelectionModel$X_nons, SelectionModel$X_rand) # joint model matrix
  mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)
  se <- sqrt(var)
  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)
  # confidence interval based on the normal approximation
  confidence_interval <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                              upper_bound = mu_hat + z * se))))
  output <- data.frame(t(data.frame("result" = c(mean = mu_hat, SE = se))))
  parameters <- matrix(c(theta_hat, theta_standard_errors),
                       ncol = 2,
                       dimnames = list(names(theta_hat),
                                       c("Estimate", "Std. Error")))
  weights_summary <- summary(as.vector(weights_nons))
  prop_scores <- c(ps_nons, est_ps_rand)


  structure(
    list(X = X,
         prop_scores = prop_scores,
         weights = as.vector(weights_nons),
         control = list(control_selection = control_selection,
                        control_outcome = control_outcome,
                        control_inference = control_inference),
         output = output,
         SE_values = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         beta = beta_statistics,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         log_likelihood = log_likelihood,
         df_residual = df_residual,
         outcome = model_out,
         selection = model_sel
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

  mu_hat <- 1/N_nons * sum(weights * weights_nons*(y - y_nons)) + 1/N_rand * sum(weights_rand * y_rand)
  mu_hat
}

