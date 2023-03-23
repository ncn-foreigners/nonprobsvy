#' nonprobIPW
#
#' nonprobIPW: Function for inference based on nonprobability big data sample and estimated propensity scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param target - `formula` with target variables.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals - an optional `named vector` with population totals.
#' @param pop_means - an optional `named vector` with population means.
#' @param pop_size - an optional `double` with population size.
#' @param method_selection - a `character` with method for propensity scores estimation
#' @param family_selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - an optional `vector` specifying strata.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
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
  weights <- rep.int(1, nrow(data)) # to remove

  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))

  # formula for outcome variable if outcome defined
  # dependents <- paste(selection, collapse = " ")
  # outcome <- stats::as.formula(paste(outcome[2], dependents))

  if (is.null(pop_totals) && !is.null(svydesign)) {
    model <- model_frame(formula = outcome, data = data, svydesign = svydesign)
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

    # Estimation for selection model
    model_sel <- internal_selection(X,
                                    X_nons,
                                    X_rand,
                                    weights,
                                    weights_rand,
                                    R,
                                    method_selection,
                                    optim_method,
                                    varcov = TRUE)

    maxLik_nons_obj <- model_sel$maxLik_nons_obj
    maxLik_rand_obj <- model_sel$maxLik_rand_obj
    log_likelihood <- model_sel$log_likelihood # maximum of the loglikelihood function
    theta_hat <- model_sel$theta
    var_cov1 <- model_sel$var_cov1
    var_cov2 <- model_sel$var_cov2

    ps_nons <- maxLik_nons_obj$ps
    est_ps_rand <- maxLik_rand_obj$ps
    hess <- maxLik_nons_obj$hess
    names(theta_hat) <- c("(Intercept)", nons_names)

    if (method_selection == "probit") { # for probit model, propensity score derivative is required
      ps_nons_der <- maxLik_nons_obj$psd
      est_ps_rand_der <- maxLik_rand_obj$psd
    }

    weights_nons <- 1/ps_nons
    N_est_nons <- sum(weights_nons)

    theta_h <- theta_h_estimation(R = R,
                                  X = X,
                                  weights_rand = weights_rand,
                                  weights = weights,
                                  h = h,
                                  method_selection = method_selection,
                                  maxit = maxit)
    names(theta_h) <- c("(Intercept)", nons_names)


    if (!is.null(pop_size)) {
      N_est_nons <- pop_size
    }

    mu_hat <- mu_hatIPW(y = y_nons,
                        weights = weights_nons,
                        N = N_est_nons) # IPW estimator

    if (control_inference$var_method == "analytic") {
      if (is.null(pop_size)) {
        b <- switch(method_selection,
                    "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - mu_hat)) %*% X_nons %*% solve(hess),
                    "cloglog" = (((1 - ps_nons)/ps_nons^2) * log(1 - ps_nons) * (y_nons - mu_hat)) %*% X_nons %*% solve(hess),
                    "probit" = - (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) %*% X_nons %*% solve(hess)
        )
      } else {
        b <- switch(method_selection,
                    "logit" = (((1 - ps_nons)/ps_nons) * y_nons) %*% X_nons %*% solve(hess),
                    "cloglog" = (((1 - ps_nons)/ps_nons^2) * log(1 - ps_nons) * y_nons) %*% X_nons %*% solve(hess),
                    "probit" = - (ps_nons_der/ps_nons^2 * (y_nons - mu_hat + 1)) %*% X_nons %*% solve(hess)
        )
      }


      # sparse matrix
      b_vec <- cbind(-1, b)
      H_mx <- cbind(0, N_est_nons * solve(hess))
      sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

      if (method_selection == "probit") {

        V1 <- var_cov1(X_nons, y_nons, mu_hat, ps_nons, ps_nons_der, pop_size) # fixed
        V2 <- var_cov2(X_rand, est_ps_rand, ps_rand, est_ps_rand_der, n_rand, N_est_nons)

      } else {

        V1 <- var_cov1(X_nons, y_nons, mu_hat, ps_nons, pop_size) # fixed
        V2 <- var_cov2(X_rand, est_ps_rand, ps_rand, n_rand, N_est_nons)

      }

      # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
      V_mx_nonprob <- sparse_mx %*% V1 %*% t(as.matrix(sparse_mx)) # nonprobability component
      V_mx_prob <- sparse_mx %*% V2 %*% t(as.matrix(sparse_mx)) # probability component - strange results for probit model
      V_mx <- V_mx_nonprob + V_mx_prob

      var_nonprob <- as.vector(V_mx_nonprob[1,1])
      var_prob <- as.vector(V_mx_prob[1,1])
      var <- as.vector(V_mx[1,1])

      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)

      # vector of variances for theta_hat
      theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))
    } else if (control_inference$var_method == "bootstrap") {
      var <- bootIPW(X_rand = X_rand,
                     X_nons = X_nons,
                     y = y_nons,
                     family_outcome = family_outcome,
                     num_boot = 1000,
                     weights = weights,
                     weights_rand = weights_rand,
                     R = R,
                     mu_hat = mu_hat,
                     method_selection = method_selection,
                     n_nons = n_nons,
                     n_rand = n_rand,
                     optim_method = optim_method,
                     pop_size = pop_size,
                     varcov = FALSE
      )
      inf <- "not computed for bootstrap variance"
    } else {
      stop("Invalid method for variance estimation.")
    }

  } else if (is.null(pop_totals) && !is.null(svydesign)) {
    # model for outcome formula
    model <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)

    theta_h <- theta_h_estimation(R = rep(1, nrow(model$X_nons)),
                                  X = model$X_nons,
                                  weights_rand = NULL,
                                  weights = weights,
                                  h = h,
                                  method_selection = method_selection,
                                  maxit = maxit,
                                  pop_totals = model$pop_totals)
    names(theta_h) <- c("(Intercept)", model$nons_names)
    method <- get_method(method_selection)
    inv_link <- method$make_link_inv
    ps_nons <- inv_link(theta_h %*% t(model$X_nons))
    N_nons <- sum(1/ps_nons)

    mu_hat <- mu_hatIPW(model$y_nons, weights = 1/ps_nons, N = N_nons)
    var <- 0
    se_nonprob <- 0
    se_prob <- 0
    theta_hat <- NULL
  } else if (is.null(pop_totals) && is.null(svydesign) && !is.null(pop_means)) {
    if (!is.null(pop_size)) pop_totals <- pop_size * pop_means

    # model for outcome formula
    model <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)

    theta_h <- theta_h_estimation(R = rep(1, nrow(model$X_nons)),
                                  X = model$X_nons,
                                  weights_rand = NULL,
                                  weights = weights,
                                  h = h,
                                  method_selection = method_selection,
                                  maxit = maxit,
                                  pop_totals = model$pop_totals)
    names(theta_h) <- c("(Intercept)", model$nons_names)
    method <- get_method(method_selection)
    inv_link <- method$make_link_inv
    ps_nons <- inv_link(theta_h %*% t(model$X_nons))
    N_nons <- sum(1/ps_nons)

    mu_hat <- mu_hatIPW(model$y_nons, weights = 1/ps_nons, N = N_nons)
    var <- 0
    se_nonprob <- 0
    se_prob <- 0
    theta_hat <- NULL
  }

  else {
    stop("Please, provide ...")
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
         SE_nonprob = ifelse(control_inference$var_method == "analytic", se_nonprob, inf),
         SE_prob = ifelse(control_inference$var_method == "analytic", se_prob, inf),
         #variance_covariance = V_mx,
         CI = ci,
         theta_h = theta_h,
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

  mu_hat <- (1/N) * sum(y * weights)
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


