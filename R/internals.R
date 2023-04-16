# These functions are only used internally in the package, so there is no need for documenting them.

internal_selection <- function(X,
                               X_nons,
                               X_rand,
                               weights,
                               weights_rand,
                               R,
                               method_selection,
                               optim_method,
                               h = h,
                               smooth,
                               maxit,
                               varcov = FALSE,
                               ...) {

  method <- get_method(method_selection)
  ps_method <- method$make_propen_score # function for propensity score estimation
  loglike <- method$make_log_like
  gradient <- method$make_gradient
  hessian <- method$make_hessian

  if (!smooth) {
    # initial values for propensity score estimation
    start <- start_fit(X,
                       R,
                       weights,
                       weights_rand,
                       method_selection)

    log_like <- loglike(X_nons,
                        X_rand,
                        weights_rand)

    gradient <- gradient(X_nons,
                         X_rand,
                         weights_rand)

    hessian <- hessian(X_nons,
                       X_rand,
                       weights_rand)

    maxLik_nons_obj <- ps_method(X_nons,
                                 log_like,
                                 gradient,
                                 hessian,
                                 start,
                                 optim_method)

    maxLik_rand_obj <- ps_method(X_rand,
                                 log_like,
                                 gradient,
                                 hessian,
                                 start,
                                 optim_method)

    theta <- maxLik_nons_obj$theta_hat
    log_likelihood <- log_like(theta)

    list(maxLik_rand_obj = maxLik_rand_obj,
         maxLik_nons_obj = maxLik_nons_obj,
         theta = theta,
         log_likelihood = log_likelihood,
         var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
         var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"))
  } else {

    inv_link <- method$make_link_inv
    h_object <- theta_h_estimation(R = R,
                                   X = X,
                                   weights_rand = weights_rand,
                                   weights = weights,
                                   h = h,
                                   method_selection = method_selection,
                                   maxit = maxit) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation
    theta_hat <- h_object$theta_h
    hess <- h_object$hess
    grad <- h_object$grad
    ps_nons <- inv_link(theta_hat %*% t(as.matrix(X_nons)))
    est_ps_rand <- inv_link(theta_hat %*% t(as.matrix(X_rand)))

    if (method_selection == "probit") { # for probit model, propensity score derivative is required
      dinv_link <- method$make_link_inv_der
      ps_nons_der <- dinv_link(theta_hat %*% t(as.matrix(X_nons)))
      est_ps_rand_der <- dinv_link(theta_hat %*% t(as.matrix(X_rand)))
    }

    list(theta_hat = theta_hat,
         hess = hess,
         grad = grad,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ifelse(method_selection == "probit", ps_nons_der, NA),
         est_ps_rand_der = ifelse(method_selection == "probit", est_ps_rand_der, NA),
         var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
         var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"))
  }

}

internal_outcome <- function(X_nons,
                             X_rand,
                             y,
                             weights,
                             family_outcome,
                             pop_totals = FALSE) {

  # estimation
  model_nons <- nonprobMI_fit(x = X_nons,
                              y = y,
                              weights = weights,
                              family_outcome = family_outcome)


  model_nons_coefs <- model_nons$coefficients

  if (pop_totals) {
    y_rand_pred <- sum(X_rand * model_nons_coefs)
  } else {
     y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample
  }
  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  list(y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       model_nons_coefs = model_nons_coefs)

}

theta_h_estimation <- function(R,
                               X,
                               weights_rand,
                               weights,
                               h,
                               method_selection,
                               maxit,
                               pop_totals = NULL,
                               pop_means = NULL){

  p <- ncol(X)
  start0 <- rep(0, p)
  # theta estimation by unbiased estimating function depending on the h_x function TODO
  u_theta <- u_theta(R = R,
                     X = X,
                     weights = c(weights_rand, weights),
                     h = h,
                     method_selection = method_selection,
                     pop_totals = pop_totals)

  u_theta_der <- u_theta_der(R = R,
                             X = X,
                             weights = c(weights_rand, weights),
                             h = h,
                             method_selection = method_selection,
                             pop_totals = pop_totals)
  for (i in 1:maxit) {
    start <- start0 + MASS::ginv(u_theta_der(start0)) %*% u_theta(start0)
    if (sum(abs(start - start0)) < 0.001) break;
    if (sum(abs(start - start0)) > 1000) break;
    start0 <- start
  }
  theta_h <- as.vector(start)

  list(theta_h = theta_h,
       hess = u_theta_der,
       grad = u_theta)
}


internal_varIPW <- function(X_nons,
                            X_rand,
                            y_nons,
                            ps_nons,
                            mu_hat,
                            hess,
                            ps_nons_der,
                            N,
                            est_ps_rand,
                            ps_rand,
                            est_ps_rand_der,
                            n_rand,
                            pop_size,
                            method_selection,
                            var_cov1 = var_cov1,
                            var_cov2 = var_cov2) {

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
  H_mx <- cbind(0, N * solve(hess))
  sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

  if (method_selection == "probit") {

    V1 <- var_cov1(X = X_nons,
                   y = y_nons,
                   mu = mu_hat,
                   ps = ps_nons,
                   psd = ps_nons_der,
                   N = pop_size) # fixed
    V2 <- var_cov2(X = X_rand,
                   eps = est_ps_rand,
                   ps = ps_rand,
                   psd = est_ps_rand_der,
                   n = n_rand,
                   N = N)

  } else {

    V1 <- var_cov1(X = X_nons,
                   y = y_nons,
                   mu = mu_hat,
                   ps = ps_nons,
                   N = pop_size) # fixed
    V2 <- var_cov2(X = X_rand,
                   eps = est_ps_rand,
                   ps = ps_rand,
                   n = n_rand,
                   N = N)

  }

  # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
  V_mx_nonprob <- sparse_mx %*% V1 %*% t(as.matrix(sparse_mx)) # nonprobability component
  V_mx_prob <- sparse_mx %*% V2 %*% t(as.matrix(sparse_mx)) # probability component
  V_mx <- V_mx_nonprob + V_mx_prob

  var_nonprob <- as.vector(V_mx_nonprob[1,1])
  var_prob <- as.vector(V_mx_prob[1,1])
  var <- as.vector(V_mx[1,1])
  # vector of variances for theta_hat
  theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

  list(var_nonprob = var_nonprob,
       var_prob = var_prob,
       var = var,
       theta_hat_var = theta_hat_var)
}


internal_varDR <- function(OutcomeModel,
                           SelectionModel,
                           y_nons_pred,
                           method_selection,
                           ps_nons,
                           hess,
                           ps_nons_der,
                           est_ps_rand,
                           y_rand_pred,
                           N_nons,
                           est_ps_rand_der,
                           svydesign) {

  h_n <- 1/N_nons * sum(OutcomeModel$y_nons - y_nons_pred) # errors mean
  b <- switch(method_selection,
              "logit" = (((1 - ps_nons)/ps_nons) * (OutcomeModel$y_nons - y_nons_pred - h_n)) %*% SelectionModel$X_nons %*% solve(hess),
              "cloglog" = (((1 - ps_nons)/ps_nons^2) * (OutcomeModel$y_nons - y_nons_pred - h_n) * log(1 - ps_nons)) %*% SelectionModel$X_nons %*% solve(hess),
              "probit" = - (ps_nons_der/ps_nons^2 * (OutcomeModel$y_nons - y_nons_pred - h_n)) %*% SelectionModel$X_nons %*% solve(hess)
  )
  # design based variance estimation based on approximations of the second-order inclusion probabilities
  t <- switch(method_selection,
              "logit" = as.vector(est_ps_rand) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_nons * sum(y_nons_pred),
              "cloglog" = as.vector(log(1 - est_ps_rand)) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_nons * sum(y_nons_pred),
              "probit" = as.vector(est_ps_rand_der/(1 - est_ps_rand)) * SelectionModel$X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_nons * sum(y_nons_pred)
  )
  svydesign <- stats::update(svydesign,
                             t = t)
  # asymptotic variance by each propensity score method (nonprobability component)
  V <- switch(method_selection,
              "logit" = (1/N_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(SelectionModel$X_nons))^2),
              "cloglog" = (1/N_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(log((1 - ps_nons)/ps_nons) * as.data.frame(SelectionModel$X_nons))))^2),
              "probit" = (1/N_nons^2) * sum((1 - ps_nons) * (((OutcomeModel$y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(ps_nons_der/(ps_nons*(1 - ps_nons)) * as.data.frame(SelectionModel$X_nons))))^2)
  )

  svydesign_mean <- survey::svymean(~t, svydesign) #perhaps using survey package to compute prob variance
  var_prob <- as.vector(attr(svydesign_mean, "var"))
  var_nonprob <- as.vector(V) #nonprob

  list(var_prob = var_prob,
       var_nonprob = var_nonprob)
}


model_frame <- function(formula, data, svydesign = NULL, pop_totals = NULL, pop_size = NULL) {

  if (!is.null(svydesign)) {
  XY_nons <- model.frame(formula, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample with intercept
  nons_names <- attr(terms(formula, data = data), "term.labels")
  if (all(nons_names %in% colnames(svydesign$variables))) {
    X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept
  } else {
    stop("variable names in data and svydesign do not match")
  }
  y_nons <- XY_nons[,1]

  list(X_nons = X_nons,
       X_rand = X_rand,
       nons_names = nons_names,
       y_nons = y_nons)
  } else if (!is.null(pop_totals)) {
    XY_nons <- model.frame(formula, data)
    X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample with intercept
    nons_names <- attr(terms(formula, data = data), "term.labels")
    if(all(nons_names %in% names(pop_totals))) { # pop_totals, pop_means defined such as in `calibrate` function
      pop_totals <-  pop_totals[nons_names]
    } else {
      warning("Selection and population totals have different names.")
    }
    y_nons <- XY_nons[,1]

    list(X_nons = X_nons,
         pop_totals = pop_totals,
         nons_names = nons_names,
         y_nons = y_nons)
  }
}


get_method <- function(method_selection) {
  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }
}
