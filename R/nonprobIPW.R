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

  weights <- rep.int(1, nrow(data)) # to remove

  # formula for outcome variable if target defined
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))

  # formula for outcome variable if outcome defined
  # dependents <- paste(selection, collapse = " ")
  # outcome <- stats::as.formula(paste(outcome[2], dependents))

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample
  X_rand <- model.matrix(selection, svydesign$variables) #matrix for probability sample
  y_nons <- XY_nons[,1]


  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand)

  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand

  # initial values for propensity score estimation
  start <- start_fit(X,
                     R,
                     weights,
                     weights_rand,
                     method_selection)


  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  ps_method <- method$make_propen_score # function for propensity score estimation
  loglike <- method$make_log_like
  gradient <- method$make_gradient
  hessian <- method$make_hessian
  var_cov1 <- method$variance_covariance1
  var_cov2 <- method$variance_covariance2

  optim_method <- control_selection$optim_method


  nonprobIPW_inference <- function(...){

    #loglike, gradient, hessian here
    log_like <- loglike(X_nons, X_rand, weights_rand)
    gradient <- gradient(X_nons, X_rand, weights_rand)
    hessian <- hessian(X_nons, X_rand, weights_rand)


    theta_hat <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$theta_hat
    hess <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$hess

    # to complete
    if(method_selection == "probit"){ # for probit model propensity score derivative is required

      ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$psd
      ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$ps
      est_ps_rand_der <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)$psd
      est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)$ps

    } else {

      ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$ps
      est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)$ps

    }


    # pearson residuals for propensity score model
    pearson_residuals <- pearson_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand)

    # deviance residuals for propensity score model
    deviance_residuals <- deviance_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand)

    weights_nons <- 1/ps_nons
    N_est_nons <- sum(weights_nons)
    # weights <- d_nons


    if (!is.null(pop_size)) {

      N_est_nons <- pop_size

    }

    mu_hat <- mu_hatIPW(y = y_nons,
                        weights = weights_nons,
                        N = N_est_nons) # IPW estimator
  #  if(method_selection == "logit") {

   #   vv <- (1-ps_nons)/ps_nons
  #    print(summary(as.vector(vv)))

   # } else if (method_selection == "cloglog") {


    #  vv <- (1 - ps_nons)/ps_nons^2 * log(1 - ps_nons)
     # print(summary(as.vector(vv)))

    #} else if (method_selection == "probit") {

     # vv <- ps_nons_der/ps_nons^2
      #print(summary(as.vector(vv)))

    #}

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
    se <- sqrt(var)

    # maximum of the  loglikelihood function
    log_likelihood <- log_like(theta_hat)

    # vector of variances for theta_hat
    theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

    # variance for mu_hat
    # var <- switch(method.selection,
    #    "logit" =  (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - mu_hat)/ps_nons) - b %*% t(as.matrix(X_nons)))^2) + D.var,
    #   "cloglog" = (1/N_est_nons^2) * sum(((1 - ps_nons)/ps_nons^2)*((y_nons - mu_hat) + b %*% t(as.matrix(X_nons)) * log(1 - ps_nons))^2) + D.var,
    #  "probit" = (1/N_est_nons^2) * sum((y_nons - mu_hat)^2 * (1 - ps_nons)/ps_nons^2 + 2 * (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) * (b %*% t(as.matrix(X_nons))) + (psdA)/((ps_nons^2)*(1 - ps_nons)) * (b %*% t(as.matrix(X_nons)))^2) + D.var,
    # )


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

      se <- sqrt(var)


    }

    mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)

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
           SE_nonprob = se_nonprob,
           SE_prob = se_prob,
           #variance_covariance = V_mx,
           CI = ci
           #theta = theta_hat,
           #theta_variance = theta_hat_var,
           #pearson_residuals = pearson_residuals,
           #deviance_residuals = deviance_residuals,
           #log_likelihood = log_likelihood
    ),
    class = "Inverse probability weighted")

  }


  infer <- nonprobIPW_inference()

  infer


}


mu_hatIPW <- function(y,
                      weights,
                      N) {

  mu_hat <- (1/N) * sum(y * weights)

  mu_hat

}

