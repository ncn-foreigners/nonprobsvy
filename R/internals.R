# These functions are only used internally in the package, so there is no need for documenting them.
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats delete.response
#' @importFrom stats model.response
#' @importFrom stats summary.glm
#' @importFrom stats contrasts
#' @importFrom nleqslv nleqslv
#' @importFrom stats get_all_vars

# Selection model object
internal_selection <- function(X,
                               X_nons,
                               X_rand,
                               weights,
                               weights_rand,
                               R,
                               method_selection,
                               optim_method,
                               h,
                               est_method,
                               maxit,
                               control_selection,
                               start,
                               bias_correction = FALSE,
                               varcov = FALSE,
                               ...) {

  if (bias_correction == TRUE) est_method <- "mm"
  estimation_method <- get_method(est_method)
  estimation_method$model_selection(X = X,
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
                                    varcov = varcov,
                                    control_selection = control_selection,
                                    start = start,
                                    ...)

}
# Outcome model object
internal_outcome <- function(outcome,
                             data,
                             weights,
                             family_outcome,
                             start_outcome) {
  # estimation
  model_nons <- nonprobMI_fit(outcome = outcome,
                              data = data,
                              weights = weights,
                              family_outcome = family_outcome,
                              start = start_outcome)
  model_nons_summary <- summary(model_nons)

  list(glm = model_nons,
       glm_summary = model_nons_summary)

}
theta_h_estimation <- function(R,
                               X,
                               weights_rand,
                               weights,
                               h,
                               method_selection,
                               maxit,
                               start = NULL,
                               pop_totals = NULL,
                               pop_means = NULL){ # TODO with BERENZ recommendation

  p <- ncol(X)
  # if (is.null(pop_totals) & is.null(pop_means)) {
  #   if (is.null(start)) {
  #     start0 <- start_fit(X = X, # <--- does not work with pop_totals
  #                         R = R,
  #                         weights = weights,
  #                         weights_rand = weights_rand,
  #                         method_selection = method_selection)
  #   } else {
  #     start0 <- start
  #   }
  # } else { # TODO customize start point for fitting with population totals
  #   # start0 <- rep(.8, ncol(X))
  #   # X_pop <- rbind(X, pop_totals)
  #   # weights_randd <- 1
  #   if (is.null(start)) {
  #     start0 <- start_fit(X = X, # <--- does not work with pop_totals
  #                         R = R,
  #                         weights = weights,
  #                         weights_rand = weights_rand,
  #                         method_selection = method_selection)
  #   } else {
  #     start0 <- start
  #   }
  # }
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

  if (method_selection == "cloglog") {
    root <- nleqslv::nleqslv(x = start,
                             fn = u_theta,
                             method = "Newton", # TODO consider the methods
                             global = "cline", #qline",
                             xscalm = "fixed",
                             jacobian = TRUE
                             )
  } else {
    root <- nleqslv::nleqslv(x = start,
                             fn = u_theta,
                             method = "Newton", # TODO consider the methods
                             global = "cline", #qline",
                             xscalm = "fixed",
                             jacobian = TRUE,
                             jac = u_theta_der
                             #control = list(sigma = 0.1, trace = 1)
    )
  }


  theta_root <- root$x
  if (root$termcd %in% c(2:7, -10)) {
    switch(as.character(root$termcd),
           "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv, but user must check if function values are acceptably small."),
           "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable new point when fitting selection model by nleqslv."),
           "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
           "5" = warning("ill-conditioned Jacobian when fitting selection model by nleqslv."),
           "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
           "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
           "-10" = warning("user specified Jacobian is incorrect when fitting selection model by nleqslv."))
  }
  theta_h <- as.vector(theta_root)
  grad <- u_theta(theta_h)
  if (method_selection == "cloglog") {
    hess <- root$jac
  } else {
    hess <- u_theta_der(theta_h) # TODO compare with root$jac
  }

  list(theta_h = theta_h,
       hess = hess,
       grad = grad)
}
# code for the function comes from the ncvreg package
setup_lambda <- function(X,
                         y,
                         weights,
                         method_selection,
                         lambda_min,
                         nlambda,
                         pop_totals,
                         alpha = 1,
                         log_lambda = FALSE,
                         ...) { #consider penalty factor here # TO consider for pop_totals/pop_means

  #fit <- glm.fit(x = X, y = y, weights = weights, family = binomial(link = method_selection))
  if (is.null(pop_totals)) {
    fit <- stats::glm(y~1,
                      weights = weights,
                      family = binomial(link = method_selection))

    n <- length(y)
    p <- ncol(X)
    w <- fit$weights
    r <- as.matrix(stats::residuals(fit, "working") * w)
    zmax <- max(crossprod(X, r))/n
    lambda_max <- zmax/alpha
  } else {
    lambda_max <- .1
  }
  if (log_lambda) { # lambda sequence on log-scale
    if (lambda_min==0) {
      lambda <- c(exp(seq(log(lambda_max), log(.001*lambda_max), length=nlambda-1)), 0)
    } else {
      lambda <- exp(seq(log(lambda_max), log(lambda_min*lambda_max), length=nlambda))
    }
  } else { # lambda sequence on linear-scale
    if (lambda_min==0) {
      lambda <- c(seq(lambda_max, 0.001*lambda_max, length = nlambda-1), 0)
    } else {
      lambda <- seq(lambda_max, lambda_min*lambda_max, length = nlambda)
    }
  }
  lambda
}

# score equation for theta, used in variable selection
u_theta <- function(R,
                    X,
                    weights,
                    method_selection,
                    h,
                    N = NULL,
                    pop_totals = NULL,
                    pop_size = NULL
) {


  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
  inv_link <- method$make_link_inv
  function(par) {
    theta <- as.matrix(par)
    n <- length(R)
    X0 <- as.matrix(X)
    eta_pi <- X0 %*% theta
    ps <- inv_link(eta_pi)
    R_rand <- 1 - R
    ps <- as.vector(ps)
    N_nons <- sum(1/ps)
    weights_sum <- sum(weights)

    if (is.null(pop_totals)) {
      eq <- switch(h,
                   "1" = c(apply(X0 * R/ps * weights - X0 * R_rand * weights, 2, sum)), # consider division by N_nons
                   "2" = c(apply(X0 * R * weights - X0 * R_rand * ps * weights, 2, sum)))
    } else {
      eq <- c(apply(X0 * R/ps * weights, 2, sum)) - pop_totals
    }
    eq
  }
}

# derivative of score equation for theta, used in variable selection
u_theta_der <-  function(R,
                         X,
                         weights,
                         method_selection,
                         h,
                         N = NULL,
                         pop_totals = NULL
)
{
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
  inv_link <- method$make_link_inv
  dinv_link <- method$make_link_inv_der
  inv_link_rev <- method$make_link_inv_rev

  function(par) {
    theta <- as.matrix(par)
    X0 <- as.matrix(X)
    p <- ncol(X0)
    eta <- X0 %*% theta
    ps <- inv_link(eta)
    ps <- as.vector(ps)
    N_nons <- sum(1/ps)
    R_rand <- 1 - R
    weights_sum <- sum(weights)

    if (!is.null(pop_totals)) {
      mxDer <- t(R * as.data.frame(X0) * weights * inv_link_rev(eta)) %*% X0
    } else {
      mxDer <-switch(h,
                     "1" = t(R * as.data.frame(X0) * weights * inv_link_rev(eta)) %*% X0, # TODO bug here when solve for some data - probably because of inv_link_rev
                     "2" = - t(R_rand * as.data.frame(X0) * weights * dinv_link(eta)) %*% X0)
    }
    as.matrix(mxDer, nrow = p) # consider division by N_nons
  }
}
# Variance for inverse probability weighted estimator
internal_varIPW <- function(svydesign,
                            X_nons,
                            X_rand,
                            y_nons,
                            weights,
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
                            pop_totals,
                            method_selection,
                            est_method,
                            theta,
                            h,
                            var_cov1 = var_cov1,
                            var_cov2 = var_cov2) {

  eta <- as.vector(X_nons %*% as.matrix(theta))
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
  b_obj <- method$b_vec_ipw(X = X_nons,
                            ps = ps_nons,
                            psd = ps_nons_der,
                            y = y_nons,
                            mu = mu_hat,
                            hess = hess,
                            eta = eta,
                            pop_size = pop_size,
                            weights = weights)
  b <- b_obj$b

  # sparse matrix
  b_vec <- cbind(-1, b)
  H_mx <- cbind(0, N * solve(hess))
  sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

  V1 <- var_cov1(X = X_nons,
                 y = y_nons,
                 mu = mu_hat,
                 ps = ps_nons,
                 psd = ps_nons_der,
                 pop_size = pop_size,
                 est_method = est_method,
                 h = h,
                 weights = weights,
                 pop_totals = pop_totals) # fixed
  V2 <- var_cov2(X = X_rand,
                 svydesign = svydesign,
                 eps = est_ps_rand,
                 est_method = est_method,
                 h = h,
                 pop_totals = pop_totals,
                 psd = est_ps_rand_der)


  # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
  V_mx_nonprob <- sparse_mx %*% V1 %*% t(as.matrix(sparse_mx)) # nonprobability component
  V_mx_prob <- sparse_mx %*% V2 %*% t(as.matrix(sparse_mx)) # probability component
  V_mx <- V_mx_nonprob + V_mx_prob

  var_nonprob <- as.vector(V_mx_nonprob[1,1])
  var_prob <- as.vector(V_mx_prob[1,1])
  var <- as.vector(V_mx[1,1])
  # vector of variances for theta_hat
  #theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

  list(var_nonprob = var_nonprob,
       var_prob = var_prob,
       var = var)
}
# Variance for doubly robust estimator
# TODO add nn and pmm
internal_varDR <- function(OutcomeModel,
                           SelectionModel,
                           y_nons_pred,
                           weights,
                           weights_rand,
                           method_selection,
                           control_selection,
                           theta,
                           ps_nons,
                           hess,
                           ps_nons_der,
                           est_ps_rand,
                           y_rand_pred,
                           N_nons,
                           est_ps_rand_der,
                           svydesign,
                           est_method,
                           h,
                           pop_totals,
                           sigma,
                           bias_correction) {

  ######### mm
  if (bias_correction == TRUE) {
    infl1 <- (weights * (OutcomeModel$y_nons - y_nons_pred))^2 / ps_nons^2
    infl2 <- (weights * (OutcomeModel$y_nons - y_nons_pred))^2 / ps_nons

    # Variance estimators ####
    svydesign <- stats::update(svydesign,
                               y_rand = y_rand_pred)
    svydesign_mean <- survey::svymean(~y_rand, svydesign)

    var_prob <- as.vector(attr(svydesign_mean, "var")) # based on survey package, probability component
    var_nonprob <- (sum((infl1) - 2*infl2) + sum(weights_rand * sigma))/N_nons^2 # TODO potential bug here nonprobability component
    } else {
    eta <- as.vector(SelectionModel$X_nons %*% as.matrix(theta))
    h_n <- 1/N_nons * sum(OutcomeModel$y_nons - y_nons_pred) # TODO add weights # errors mean
    method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
    method <- get_method(method_selection)
    est_method <- get_method(est_method)
    #psd <- method$make_link_inv_der(eta)

    b <- method$b_vec_dr(X = SelectionModel$X_nons,
                         ps = ps_nons,
                         psd = ps_nons_der,
                         y = OutcomeModel$y_nons,
                         hess = hess,
                         eta = eta,
                         h_n = h_n,
                         y_pred = y_nons_pred,
                         weights = weights)

    # asymptotic variance by each propensity score method (nonprobability component)
    var_nonprob <- est_method$make_var_nonprob(ps = ps_nons,
                                               psd = ps_nons_der,
                                               y = OutcomeModel$y_nons,
                                               y_pred = y_nons_pred,
                                               h_n = h_n,
                                               X = SelectionModel$X_nons,
                                               b = b,
                                               N = N_nons,
                                               h = h,
                                               method_selection = method_selection,
                                               weights = weights,
                                               pop_totals = pop_totals)


    if (is.null(pop_totals)) {
      t <- est_method$make_t(X = SelectionModel$X_rand,
                             ps = est_ps_rand,
                             psd = est_ps_rand_der,
                             b = b,
                             h = h,
                             y_rand = y_rand_pred,
                             y_nons = y_nons_pred,
                             N = N_nons,
                             method_selection = method_selection,
                             weights = weights)
      # design based variance estimation based on approximations of the second-order inclusion probabilities
      svydesign <- stats::update(svydesign,
                                 t = t)
      svydesign_mean <- survey::svymean(~t, svydesign) #perhaps using survey package to compute prob variance
      var_prob <- as.vector(attr(svydesign_mean, "var"))
    } else {
      var_prob <- 0
    }
  }

  list(var_prob = var_prob,
       var_nonprob = var_nonprob)
}
# Variance for mass imputation estimator
internal_varMI <- function(svydesign,
                           X_nons,
                           X_rand,
                           y,
                           y_pred,
                           y_hat,
                           weights_rand,
                           method,
                           n_rand,
                           n_nons,
                           N,
                           family,
                           parameters,
                           pop_totals
                           ) {

  if(is.character(family)) {
    family_nonprobsvy <- paste(family, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (is.null(pop_totals)) {
    svydesign_mean <- survey::svymean(~y_hat_MI, svydesign)
    var_prob <- as.vector(attr(svydesign_mean, "var")) # probability component, should be bigger for nn
    if (method == "nn") {
      sigma_hat <- mean((y - y_pred)^2) # family_nonprobsvy$variance(mu = y_pred, y  = y)
      est_ps  <- n_nons/N
      var_nonprob <- n_rand/N^2 * (1 - est_ps)/est_ps * sigma_hat
    } else if (method == "glm") { # TODO add variance for count binary outcome variable control_outcome$method

      beta <- parameters[,1]
      eta_nons <- X_nons %*% beta
      eta_rand <- X_rand %*% beta

      mx <- 1/N * colSums(as.data.frame(X_rand) * (weights_rand * family_nonprobsvy$mu_der(eta_rand)))
      c <- solve(1/n_nons * t(as.data.frame(X_nons) * family_nonprobsvy$mu_der(eta_nons)) %*% X_nons) %*% mx
      residuals <- family_nonprobsvy$residuals(mu = y_pred, y  = y)

      # nonprobability component
      var_nonprob <- 1/n_nons^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2
      var_nonprob <- as.vector(var_nonprob)
    } else if (method == "pmm") {

      # beta <- parameters[,1]
      # eta_nons <- X_nons %*% beta
      # eta_rand <- X_rand %*% beta
      #
      # mx <- 1/N * colSums(as.data.frame(X_rand) * (weights_rand * family_nonprobsvy$mu_der(eta_rand)))
      # c <- solve(1/n_nons * t(as.data.frame(X_nons) * family_nonprobsvy$mu_der(eta_nons)) %*% X_nons) %*% mx
      # residuals <- family_nonprobsvy$residuals(mu = y_pred, y  = y)
      #
      # # nonprobability component
      # var_nonprob <- 1/n_nons^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2
      # var_nonprob <- as.vector(var_nonprob)

      # nonprobability component
      # var_nonprob <- 1/n_nons^2 * residuals^2 * X_nons %*% t(X_nons)
      var_nonprob <- 0
      # var_nonprob <- as.vector(var_nonprob)
      # TODO to consider
    }
  } else {
    if (method == "nn") {
      sigma_hat <- mean((y - y_pred)^2) # family_nonprobsvy$variance(mu = y_pred, y  = y)
      est_ps  <- n_nons/N
      var_nonprob <- n_nons/N^2 * (1 - est_ps)/est_ps * sigma_hat # what instead of n_rand here (?) now just n_nons
    } else if (method == "glm") {
      beta <- parameters[,1]
      eta_nons <- X_nons %*% beta
      if (family %in% c("binomial", "poisson")) { # TODO consider this chunk of code
        eta_rand <- pop_totals %*% beta / pop_totals[1]
      } else {
        eta_rand <- pop_totals %*% beta
      }
      mx <- 1/N * pop_totals * as.vector(family_nonprobsvy$mu_der(eta_rand))
      c <- solve(1/n_nons * t(as.data.frame(X_nons) * family_nonprobsvy$mu_der(eta_nons)) %*% X_nons) %*% mx
      residuals <- family_nonprobsvy$residuals(mu = y_pred, y  = y)

      # nonprobability component
      var_nonprob <- 1/n_nons^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2
      var_nonprob <- as.vector(var_nonprob)
    } else if (method == "pmm") {

      # beta <- parameters[,1]
      # eta_nons <- X_nons %*% beta
      #
      # if (family %in% c("binomial", "poisson")) { # TODO consider this chunk of code
      #   eta_rand <- pop_totals %*% beta / pop_totals[1]
      # } else {
      #   eta_rand <- pop_totals %*% beta
      # }
      #
      # residuals <- family_nonprobsvy$residuals(mu = y_pred, y  = y)

      # nonprobability component
      # var_nonprob <- 1/n_nons^2 * t(as.matrix(residuals^2)) %*% (family_nonprobsvy$mu_der(eta_nons) %*% t(X_nons))^2
      var_nonprob <- 0
      var_nonprob <- as.vector(var_nonprob)

    }
    var_prob <- 0
  }

  list(var_prob = var_prob,
       var_nonprob = var_nonprob)
}
# create an object with model frames and matrices to preprocess
model_frame <- function(formula, data, weights = NULL, svydesign = NULL, pop_totals = NULL, pop_size = NULL, flag = TRUE) {

  if (!is.null(svydesign)) {
  ##### Model frame for nonprobability sample #####
  model_Frame <- model.frame(formula, data)
  y_nons <- model.response(model_Frame)
  outcome_name <- names(model_Frame)[1]
  mt <- attr(model_Frame, "terms")
  nons_names <- attr(mt, "term.labels") # colnames(get_all_vars(formula, data)) names of variables of nonprobability sample terms(formula, data = data)
  ##### Model frame for probability sample #####
  if (outcome_name %in% colnames(svydesign$variables)) {
    design_to_frame <- svydesign$variables
    design_to_frame[,outcome_name][is.na(design_to_frame[,outcome_name])] <- 0
    model_Frame_rand <- model.frame(formula, design_to_frame)
    mt_rand <- attr(model_Frame_rand, "terms")
    nons_names_rand <- attr(mt_rand, "term.labels")
  } else {
    model_Frame_rand <- model.frame(formula[-2], svydesign$variables)
    mt_rand <- attr(model_Frame_rand, "terms")
    nons_names_rand <- attr(mt_rand, "term.labels")
  }
  #print(nons_names_rand)
  if (all(nons_names %in% nons_names_rand)) { #colnames(svydesign$variables)
    dot_check <- sapply(formula, FUN = function(x) {x == "."})
    if (length(formula) == 2) nons_names <- nons_names[-1]
    if (any(dot_check)) {
      xx <- paste("~", paste(nons_names, collapse = "+"))
      formula <- as.formula(paste(outcome_name, xx))
      X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables[, nons_names])
    } else {
      X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables) #matrix of probability sample with intercept
    }
    frame_nons <- model.frame(formula, data)
    X_nons <- model.matrix(frame_nons, data) #matrix for nonprobability sample with intercept
    #if (outcome) {
    #  xx <- paste("~", paste(nons_names[2:length(nons_names)], collapse = "+"))
    #  formula <- as.formula(paste(formula[2], xx))
    #  X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables[, nons_names])
    #} else {
    #  xx <- paste("~", paste(nons_names, collapse = "+"))
    #  formula <- as.formula(xx)
    #  X_rand <- model.matrix(formula, svydesign$variables[, nons_names])# matrix of probability sample with intercept
    #  }
    } else {
    stop("Variable names in data and svydesign do not match")
    }

  list(X_nons = X_nons,
       X_rand = X_rand,
       nons_names = nons_names,
       y_nons = y_nons,
       outcome_name = outcome_name,
       model_frame_rand = model_Frame_rand)

  } else if (!is.null(pop_totals)) {
    model_Frame <- model.frame(formula, data)
    X_nons <- model.matrix(model_Frame, data)
    #matrix for nonprobability sample with intercept
    #X_nons <- model.matrix(XY_nons, data, contrasts.arg = list(klasa_pr = contrasts(as.factor(XY_nons[,dep_name]), contrasts = FALSE)))
    #nons_names <- attr(terms(formula, data = data), "term.labels")
    #nons_names <- colnames(X_nons)
    #pop_totals <- pop_totals[which(attr(X_nons, "assign") == 1)]
    mt <- attr(model_Frame, "terms")
    #nons_names <- attr(mt, "term.labels")
    total_names <- colnames(X_nons)
    if (flag) {
      if(all(total_names %in% names(pop_totals))) { # TODO verify whether this warming works well.. pop_totals, pop_means defined such as in `calibrate` function
        pop_totals <- pop_totals[total_names]
      } else {
        warning("Selection and population totals have different names.")
      }
    }
    y_nons <- model.response(model_Frame)
    outcome_name <- names(model_Frame)[1]

    list(X_nons = X_nons,
         pop_totals = pop_totals,
         total_names = total_names,
         y_nons = y_nons,
         outcome_name = outcome_name,
         X_rand = NULL)
  }
}

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
                                control = list(epsilon = control_selection$epsilon,
                                               maxit = control_selection$maxit,
                                               trace = control_selection$trace)

  )
  start_model$coefficients
}

# Function for getting function from the selected method
get_method <- function(method) {
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }
  method
}

# summary helper functions
# for now just a rough sketch
specific_summary_info <- function(object, ...) {
  UseMethod("specific_summary_info")
}

specific_summary_info.nonprobsvy_ipw <- function(object,
                                                 ...) {

  coeffs_sel <- matrix(c(object$selection$coefficients, object$selection$std_err),
                  ncol = 2,
                  dimnames = list(names(object$selection$coefficients),
                                  c("Estimate", "Std. Error")))
  res <- list(
    coeffs_sel = coeffs_sel,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )

  attr(res$coeffs_sel, "glm") <- TRUE
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE # TODO
  attr(res, "model")     <- c("glm regression on selection variable")
  res
}

specific_summary_info.nonprobsvy_mi <- function(object,
                                                ...) {

  if (object$outcome[[1]]$method == "glm") { # TODO for pmm
    coeffs_out <- matrix(c(object$outcome[[1]]$coefficients, object$outcome[[1]]$std_err),
                   ncol = 2,
                   dimnames = list(names(object$outcome[[1]]$coefficients),
                                   c("Estimate", "Std. Error")))
  } else {
    coeffs_out <- "no coefficients"
  }

  res <- list(
    coeffs_out = coeffs_out
  )
  if (object$outcome[[1]]$method == "glm") {
  attr(res$coeffs_out, "glm") <- TRUE
  attr(res, "model") <- "glm regression on outcome variable"
  } else if (object$outcome[[1]]$method == "nn") {
    attr(res$coeffs_out, "glm") <- FALSE
  } else if (object$outcome[[1]]$method == "pmm") { # TODO
    attr(res$coeffs_out, "glm") <- FALSE
    # attr(res, "model") <- "glm regression on outcome variable"
  }
  res
}

specific_summary_info.nonprobsvy_dr <- function(object,
                                                ...) {

  coeffs_sel <- matrix(c(object$selection$coefficients, object$selection$std_err),
                  ncol = 2,
                  dimnames = list(names(object$selection$coefficients),
                                  c("Estimate", "Std. Error")))


  if (object$outcome[[1]]$method == "glm") {
    coeffs_out <- matrix(c(object$outcome[[1]]$coefficients, object$outcome[[1]]$std_err),
                   ncol = 2,
                   dimnames = list(names(object$outcome[[1]]$coefficients),
                                   c("Estimate", "Std. Error")))
  } else {
    coeffs_out <- "no coefficients"
  }

  res <- list(
    coeffs_sel = coeffs_sel,
    coeffs_out  = coeffs_out,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )
  attr(res$coeffs_sel, "glm") <- TRUE
  if (object$outcome[[1]]$method == "glm") {
    attr(res$coeffs_out, "glm") <- TRUE
    attr(res, "model")     <- c("glm regression on selection variable",
                                "glm regression on outcome variable")
  } else if (object$outcome[[1]]$method == "nn") {
    attr(res$coeffs_out, "glm") <- FALSE
    attr(res, "model")     <- c("glm regression on selection variable")
  }
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE

  res
}

ff <- function(formula) {
  fff <- as.character(formula)
  f <- strsplit(fff[2], "\\s*\\+\\s*")[[1]]
  outcome_formulas <- list()
  if(any(duplicated(f))) {
    warning("No unique names of the outcome variables in formula. The error has been corrected")
    f <- unique(f)
  }
  l <- length(f)
  for (i in 1:l) {
    outcome_formulas[[i]] <-  as.formula(paste(f[i], fff[3], sep = " ~ "))
  }

  list(f = f,
       outcomes = outcome_formulas,
       l = l)
}
