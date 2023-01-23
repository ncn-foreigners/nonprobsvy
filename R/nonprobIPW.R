#' nonprobIPW
#
#' nonprobIPW: Function for inference based on nonprobability big data sample and estimated propensity scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop.totals - an optional `named vector` with population totals.
#' @param pop.means - an optional `named vector` with population means.
#' @param method.selection - a `character` with method for propensity scores estimation
#' @param family.selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action a
#' @param control.selection a
#' @param control.inference a
#' @param start a
#' @param verbose a
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @export



nonprobIPW <- function(selection,
                       data,
                       svydesign,
                       pop.totals,
                       pop.means,
                       method.selection,
                       family.selection = "binomial",
                       subset,
                       weights,
                       na.action,
                       control.selection = controlSel(),
                       control.inference = controlInf(),
                       start,
                       verbose,
                       contrasts,
                       model,
                       x,
                       y,
                       ...){

  weights <- rep.int(1, nrow(data)) # to remove

  y_name <- colnames(data)[!colnames(data) %in% colnames(svydesign$variables)]
  dependents <- colnames(data)[colnames(data) %in% colnames(svydesign$variables)]
  outcome <- as.formula(paste(paste(y_name, "~"), paste(dependents, collapse = "+"))) # formula for outcome variable
  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data)
  X_rand <- model.matrix(selection, svydesign$variables)
  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand

  start <- start.fit(X_nons, #initial values for Propensity score estimation
                     X_rand,
                     weights,
                     d_rand,
                     method.selection)

  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  ps_method <- method$PropenScore
  loglike <- method$MakeLogLike
  gradient <- method$MakeGradient
  hessian <- method$MakeHessian
  VarCov1 <- method$VarianceCov1
  VarCov2 <- method$VarianceCov2

  nonprobIPW.inference <- function(...){

    #loglike, gradient, hessian here
    log_like <- loglike(X_rand, X_nons, d_rand)
    gradient <- gradient(X_rand, X_nons, d_rand)
    hessian <- hessian(X_rand, X_nons, d_rand)

    n_rand <- nrow(X_rand)

    theta_hat <- ps_method(X_nons, log_like, gradient, hessian, start)$theta_hat
    hess <- ps_method(X_nons, log_like, gradient, hessian, start)$hess

    # to complete
    if(method.selection == "probit"){

      ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian, start)$psd
      ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start)$ps
      est_ps_rand_der <- ps_method(X_rand, log_like, gradient, hessian, start)$psd
      est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start)$ps

    } else {

      ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start)$ps
      est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start)$ps

    }

    d_nons <- 1/ps_nons
    N_est_nons <- sum(d_nons)

    mu_hat <- (1/N_est_nons) * sum(y_nons/ps_nons)


    b <- switch(method.selection,
                "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - mu_hat)) %*% X_nons %*% solve(hess),
                "cloglog" = (((1 - ps_nons)/ps_nons^2) * log(1 - ps_nons) * (y_nons - mu_hat)) %*% X_nons %*% solve(hess),
                "probit" = (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) %*% X_nons %*% solve(hess)
    )


    b_vec <- cbind(-1, b)
    H_mx <- cbind(0, - N_est_nons * solve(hess))
    sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

    ##
    if(method.selection == "probit"){

      V1 <- VarCov1(X_nons, y_nons, mu_hat, ps_nons, ps_nons_der, N_est_nons) # to fix

      V2 <- VarCov2(X_rand, est_ps_rand, ps_rand, est_ps_rand_der, b, n_rand, N_est_nons)

    } else {

      V1 <- VarCov1(X_nons, y_nons, mu_hat, ps_nons, N_est_nons) # to fix

      V2 <- VarCov2(X_rand, est_ps_rand, ps_rand, b, n_rand, N_est_nons)

    }

    # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
    V_mx <- sparse_mx %*% (V1 + V2) %*% t(as.matrix(sparse_mx))


    theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

    # variance for mu_hat
    # var <- switch(method.selection,
    #    "logit" =  (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - mu_hat)/ps_nons) - b %*% t(as.matrix(X_nons)))^2) + D.var,
    #   "cloglog" = (1/N_est_nons^2) * sum(((1 - ps_nons)/ps_nons^2)*((y_nons - mu_hat) + b %*% t(as.matrix(X_nons)) * log(1 - ps_nons))^2) + D.var,
    #  "probit" = (1/N_est_nons^2) * sum((y_nons - mu_hat)^2 * (1 - ps_nons)/ps_nons^2 + 2 * (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) * (b %*% t(as.matrix(X_nons))) + (psdA)/((ps_nons^2)*(1 - ps_nons)) * (b %*% t(as.matrix(X_nons)))^2) + D.var,
    # )


    ci <- c(mu_hat - 1.96 * sqrt(V_mx[1,1]), mu_hat + 1.96 * sqrt(V_mx[1,1]))


    return(list("Population mean estimator" = mu_hat,
                "variance for mu_hat" = V_mx[1,1],
                "variance-covariance" = V_mx,
                "CI" = ci,
                "theta_hat" = theta_hat,
                "variance for theta_hat" = theta_hat_var
    ))


  }


  infer <- nonprobIPW.inference()


  return(infer)


}

