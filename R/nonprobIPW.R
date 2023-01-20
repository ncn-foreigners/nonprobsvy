#' nonprobIPW
#
#' nonprobIPW: Function for inference based on nonprobability big data sample and estimated propensity scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop.totals - an optional `named vector` with population totals.
#' @param pop.means - an optional `named vector` with population means.
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
#' @export



nonprobIPW <- function(selection,
                       data,
                       svydesign,
                       pop.totals,
                       pop.means,
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


  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data)
  X_rand <- model.matrix(outcome, svydesign$variables)
  y_nons = XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  ps_method <- method.selection$PropenScore
  loglike <- method.selection$MakeLogLike
  gradient <- method.selection$MakeGradient
  hessian <- method.selection$MakeHessian


  infer <- nonprobIPW.inference()


  return(infer)


}


nonprobIPW.inference <- function(...){

  #loglike, gradient, hessiann here

  log_like <- log_like(X_rand, X_nons, d_rand)
  gradient <- gradient(X_rand, X_nons, d_rand)
  hessian <- hessian(X_rand, X_nons, d_rand)

    ps_nons <- ps_method(X_nons, log_like, gradient, hessian)$ps
    # to complete
    if(method.selection == "probit"){

    ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian)$psd
    ps_nons <- ps_method(X_nons, log_like, gradient, hessian)$ps

    }

    d_nons <- 1/ps_nons
    N_est_nons <- sum(d_nons)

    mu_hat <- (1/N_est_nons) * sum(XY_nons[,1]/ps_nons)

    n_rand <- nrow(X_rand)

    hess <- ps_method(X_nons, log_like, gradient, hessian)$hess

    b <- switch(method.selection,
           "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - mu_hat)) %*% as.matrix(X_nons) %*% solve(hess),
           "cloglog" = (((1 - ps_nons)/ps_nons^2) * log(1 - ps_nons) * (y_nons - mu_hat)) %*% as.matrix(X_nons) %*% solve(hess),
           "probit" = (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) %*% as.matrix(X_nons) %*% solve(hess)
           )

    est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian)$ps

    s <- est_ps_rand * X_rand
    ci <- n_rand/(n_rand-1) * (1 - est_ps_rand)
    B_hat <- (t(as.matrix(ci)) %*% as.matrix(s/est_ps_rand))/sum(ci)
    ei <- (s/est_ps_rand) - B_hat
    db_var <- t(as.matrix(ei * ci)) %*% as.matrix(ei)

    D <- (1/N_est^2) * db_var
    D.var <- b %*% D %*% t(b)

    var <- switch(method.selection,
                 "logit" =  (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - mu_hat)/ps_nons) - b %*% t(as.matrix(X_nons)))^2) + D.var,
                 "cloglog" = (1/N_est_nons^2) * sum(((1 - ps_nons)/ps_nons^2)*((y_nons - mu_hat) + b %*% t(as.matrix(X_nons)) * log(1 - ps_nons))^2) + D.var,
                 "probit" = (1/N_est_nons^2) * sum((y_nons - mu_hat)^2 * (1 - ps_nons)/ps_nons^2 + 2 * (ps_nons_der/ps_nons^2 * (y_nons - mu_hat)) * (b %*% t(as.matrix(X_nons))) + (psdA)/((ps_nons^2)*(1 - ps_nons)) * (b %*% t(as.matrix(X_nons)))^2) + D.var,
                 )


    ci <- c(mu_hat - 1.96 * sqrt(var), mu_hat + 1.96 * sqrt(var))



  return(list("Population mean estimator" = mu_hat, "variance" = var, "CI" = ci))


}
