#' nonprobDR
#
#' nonprobDR: Function for inference based on nonprobability surveys, nonprobabilty big data sample and estimated propensoty scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals - an optional `named vector` with population totals.
#' @param pop_means - an optional `named vector` with population means.
#' @param pop_size - an optional `double` with population size.
#' @param method_selection - a `character` with method for propensity scores estimation
#' @param method_outcome - a `character` with method for response variable estimation
#' @param family_selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - an optional `vector` specifying strata.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
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
  weights <- rep.int(1, nrow(data)) # to remove


  if (is.null(pop_totals) && !is.null(svydesign)) {
    # model for outcome formula
    OutcomeModel <- model_frame(formula = outcome, data = data, svydesign = svydesign)

    #model for selection formula
    SelectionModel <- model_frame(formula = selection, data = data, svydesign = svydesign)
    X_sel <- rbind(SelectionModel$X_rand, SelectionModel$X_nons)
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
    model_out <- internal_outcome(OutcomeModel$X_nons,
                                  OutcomeModel$X_rand,
                                  OutcomeModel$y_nons,
                                  weights,
                                  family_outcome)

    y_rand_pred <- model_out$y_rand_pred
    y_nons_pred <- model_out$y_nons_pred
    model_nons_coefs <- model_out$model_nons_coefs

    # Estimation for selection model
    model_sel <- internal_selection(X_sel,
                                    SelectionModel$X_nons,
                                    SelectionModel$X_rand,
                                    weights,
                                    weights_rand,
                                    R,
                                    method_selection,
                                    optim_method)

    maxLik_nons_obj <- model_sel$maxLik_nons_obj
    maxLik_rand_obj <- model_sel$maxLik_rand_obj
    log_likelihood <- model_sel$log_likelihood # maximum of the loglikelihood function
    theta_hat <- model_sel$theta

    ps_nons <- maxLik_nons_obj$ps
    est_ps_rand <- maxLik_rand_obj$ps
    hess <- maxLik_nons_obj$hess
    names(theta_hat) <- c("(Intercept)", SelectionModel$nons_names)

    if (method_selection == "probit") { # for probit model, propensity score derivative is required
      ps_nons_der <- maxLik_nons_obj$psd
      est_ps_rand_der <- maxLik_rand_obj$psd
    }

    weights_nons <- 1/ps_nons
    N_est_nons <- sum(weights_nons)
    N_est_rand <- sum(weights_rand)

   # if(!is.null(pop_size)){

   #  N_est_rand <- pop_size
   #  N_est_nons <- pop_size
   # }

    mu_hat <- mu_hatDR(y = OutcomeModel$y_nons,
                       y_nons = y_nons_pred,
                       y_rand = y_rand_pred,
                       weights_nons = weights_nons,
                       weights_rand = weights_rand,
                       N_nons = N_est_nons,
                       N_rand = N_est_rand) #DR estimator

    # updating probability sample by adding y_hat variable
    svydesign <- stats::update(svydesign,
                               .y_hat_MI = y_rand_pred)


    theta_h <- theta_h_estimation(R = R,
                                  X = X_sel,
                                  weights_rand = weights_rand,
                                  weights = weights,
                                  h = h,
                                  method_selection = method_selection,
                                  maxit = maxit)
    names(theta_h) <- c("(Intercept)", SelectionModel$nons_names)

    if (control_inference$var_method == "analytic") {

      h_n <- 1/N_est_nons * sum(OutcomeModel$y_nons - y_nons_pred) # errors mean
      b <- switch(method_selection,
                  "logit" = (((1 - ps_nons)/ps_nons) * (OutcomeModel$y_nons - y_nons_pred - h_n)) %*% SelectionModel$X_nons %*% solve(hess),
                  "cloglog" = (((1 - ps_nons)/ps_nons^2) * (OutcomeModel$y_nons - y_nons_pred - h_n) * log(1 - ps_nons)) %*% SelectionModel$X_nons %*% solve(hess),
                  "probit" = - (ps_nons_der/ps_nons^2 * (OutcomeModel$y_nons - y_nons_pred - h_n)) %*% SelectionModel$X_nons %*% solve(hess)
      )
      # design based variance estimation based on approximations of the second-order inclusion probabilities
      t <- switch(method_selection,
                  "logit" = as.vector(est_ps_rand) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                  "cloglog" = as.vector(log(1 - est_ps_rand)) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                  "probit" = as.vector(est_ps_rand_der/(1 - est_ps_rand)) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred)
      )
     svydesign <- stats::update(svydesign,
                                 t = t)
      # asymptotic variance by each propensity score method (nonprobability component)
      V <- switch(method_selection,
                  "logit" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(SelectionModel$X_nons))^2),
                  "cloglog" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(log((1 - ps_nons)/ps_nons) * as.data.frame(SelectionModel$X_nons))))^2),
                  "probit" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(ps_nons_der/(ps_nons*(1 - ps_nons)) * as.data.frame(SelectionModel$X_nons))))^2)
      )

      svydesign_mean <- survey::svymean(~t, svydesign) #perhaps using survey package to compute prob variance
      var_prob <- as.vector(attr(svydesign_mean, "var"))
      var_nonprob <- as.vector(V) #nonprob

      se_prob <- sqrt(var_prob)
      se_nonprob <- sqrt(V)

      var <- var_prob + var_nonprob

    } else if (control_inference$var_method == "bootstrap") {
      var <- bootDR(SelectionModel = SelectionModel,
                    OutcomeModel = OutcomeModel,
                    family_outcome = family_outcome,
                    num_boot = 1000,
                    weights = weights,
                    weights_rand = weights_rand,
                    R = R,
                    mu_hat = mu_hat,
                    method_selection = method_selection,
                    n_nons = n_nons,
                    n_rand = n_rand,
                    optim_method = optim_method
                    )
      inf <- "not computed for bootstrap variance"
    } else {
      stop("Invalid method for variance estimation.")
    }

  } else if (is.null(pop_totals) && !is.null(svydesign)) {
    # model for outcome formula
    OutcomeModel <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)
    #model for selection formula
    SelectionModel <- model_frame(formula = selection, data = data, pop_totals = pop_totals)
    theta_h <- theta_h_estimation(R = rep(1, nrow(SelectionModel$X_nons)),
                                  X = SelectionModel$X_nons,
                                  weights_rand = NULL,
                                  weights = weights,
                                  h = h,
                                  method_selection = method_selection,
                                  maxit = maxit,
                                  pop_totals = SelectionModel$pop_totals)
    names(theta_h) <- c("(Intercept)", SelectionModel$nons_names)
    method <- get_method(method_selection)
    inv_link <- method$make_link_inv
    ps_nons <- inv_link(theta_h %*% t(SelectionModel$X_nons))
    N_est <- sum(1/ps_nons)

    model_out <- internal_outcome(X_nons = OutcomeModel$X_nons,
                                  X_rand = c(pop_size, OutcomeModel$pop_totals),
                                  y = OutcomeModel$y_nons,
                                  weights = weights,
                                  family_outcome = family_outcome)

    y_rand_pred <- model_out$y_rand_pred
    y_nons_pred <- model_out$y_nons_pred
    model_nons_coefs <- model_out$model_nons_coefs
    print(y_rand_pred)

    mu_hat <- 1/N_est * sum((1/ps_nons)*(OutcomeModel$y_nons - y_nons_pred)) + 1/pop_size * y_rand_pred

    var <- 0
    se_nonprob <- 0
    se_prob <- 0
    theta_hat <- NULL

  } else {
    stop("Please, provide ...")
  }

  se <- sqrt(var)
  alpha <- control_inference$alpha
  z <- stats::qnorm(1-alpha/2)

  # confidence interval based on the normal approximation
  ci <- c(mu_hat - z * se, mu_hat + z * se)

  # case when samples overlap - to finish
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

  mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)

  structure(
    list(mean = mu_hat,
         #VAR = var,
         #VAR_nonprob = V,
         #VAR_prob = W,
         SE = se,
         SE_nonprob = ifelse(control_inference$var_method == "analytic", se_nonprob, inf),
         SE_prob = ifelse(control_inference$var_method == "analytic", se_prob, inf),
         CI = ci,
         theta_h = theta_h,
         theta = theta_hat,
         #pearson_residuals = pearson_residuals,
         #deviance_residuals = deviance_residuals,
         #log_likelihood = log_likelihood,
         beta = model_nons_coefs
         ),
    class = "Doubly-robust")
}

#' mu_hatDR
#
#' mu_hatDR: Function for outcome variable estimation based on doubly robust estimation
#'
#' @param y - a
#' @param y_nons - a
#' @param y_rand - a
#' @param weights_nons - a
#' @param weights_rand - a
#' @param N_nons - a
#' @param N_rand - a

mu_hatDR <- function(y,
                     y_nons,
                     y_rand,
                     weights_nons,
                     weights_rand,
                     N_nons,
                     N_rand) {

  mu_hat <- 1/N_nons * sum(weights_nons*(y - y_nons)) + 1/N_rand * sum(weights_rand * y_rand)

  mu_hat

}

