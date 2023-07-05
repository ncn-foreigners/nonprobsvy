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
                               h = h,
                               est_method,
                               maxit,
                               varcov = FALSE,
                               ...) {

  estimation_method <- get_method(est_method)
  estimation_method$model_selection(X,
                                    X_nons,
                                    X_rand,
                                    weights,
                                    weights_rand,
                                    R,
                                    method_selection,
                                    optim_method,
                                    h = h,
                                    est_method,
                                    maxit,
                                    varcov,
                                    ...)

}
# Outcome model object
internal_outcome <- function(outcome,
                             data,
                             weights,
                             family_outcome) {
  # estimation
  model_nons <- nonprobMI_fit(outcome = outcome,
                              data = data,
                              weights = weights,
                              family_outcome = family_outcome)
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
                               pop_totals = NULL,
                               pop_means = NULL){

  p <- ncol(X)
  start0 <- start_fit(X = X, # <--- does not work with pop_totals
                      R = R,
                      weights = weights,
                      weights_rand = weights_rand,
                      method_selection = method_selection)
  #start0 <- rep(0, p)
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

  #root <- rootSolve::multiroot(u_theta,
  #                             jacfunc = u_theta_der,
  #                             start = start0)
  #print(root$root)
  root <- nleqslv::nleqslv(x = start0,
                           fn = u_theta,
                           method = "Newton", # TODO consider the methods
                           global = "qline",
                           xscalm = "fixed",
                           jacobian = TRUE,
                           jac = u_theta_der
                           )
  start <- root$x
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
  #it <- 0
  #for (i in 1:maxit) {
  #  it <- it + 1
  #  start <- start0 - MASS::ginv(u_theta_der(start0)) %*% u_theta(start0) # consider solve function
  #  if (sum(abs(start - start0)) < 0.001) break;
  #  if (sum(abs(start - start0)) > 1000)  {
  #    warning("algorithm did not converge")
  #    break
  #  }
  #  start0 <- start
  #}
  #if (it == maxit) warning("algorithm did not converge - iteration limit exceeded.")
  theta_h <- as.vector(start)
  grad <- u_theta(theta_h)
  hess <- u_theta_der(theta_h) # TODO compare with root$jac
  #hess <- root$jac

  list(theta_h = theta_h,
       hess = hess,
       grad = grad)
}
# Variance for inverse probability weighted estimator
internal_varIPW <- function(svydesign,
                            selection_formula,
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
  hess_inv <- b_obj$hess_inv

  # sparse matrix
  b_vec <- cbind(-1, b)
  H_mx <- cbind(0, N * hess_inv)
  sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

  V1 <- var_cov1(X = X_nons,
                 y = y_nons,
                 mu = mu_hat,
                 ps = ps_nons,
                 psd = ps_nons_der,
                 pop_size = pop_size,
                 est_method = est_method,
                 h = h,
                 weights = weights) # fixed
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
internal_varDR <- function(OutcomeModel,
                           SelectionModel,
                           y_nons_pred,
                           weights,
                           method_selection,
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
                           pop_totals) {

  eta <- as.vector(SelectionModel$X_nons %*% as.matrix(theta))
  h_n <- 1/N_nons * sum(OutcomeModel$y_nons - y_nons_pred) # TODO add weights # errors mean
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
                           family
                           ) {

  svydesign_mean <- survey::svymean(~y_hat_MI, svydesign)
  var_prob <- as.vector(attr(svydesign_mean, "var")) # probability component, should be bigger for nn

  if (method == "nn") {
    if(is.character(family)) {
      family_nonprobsvy <- paste(family, "_nonprobsvy", sep = "")
      family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
      family_nonprobsvy <- family_nonprobsvy()
    }

    sigma_hat <- family_nonprobsvy$variance(mu = y_pred, y  = y)
    est_ps  <- n_nons/N
    var_nonprob <- n_rand/N^2 * (1 - est_ps)/est_ps * sigma_hat

  } else if (method == "glm") { # control_outcome$method

    mx <- 1/N * colSums(weights_rand * X_rand)
    c <- solve(1/n_nons * t(X_nons) %*% X_nons) %*% mx
    e <- y - y_pred

    # nonprobability component
    var_nonprob <- 1/n_nons^2 * t(as.matrix(e^2)) %*% (X_nons %*% c)^2
    var_nonprob <- as.vector(var_nonprob)
  }

  list(var_prob = var_prob,
       var_nonprob = var_nonprob)
}
# create an object with model frames and matrices to preprocess
model_frame <- function(formula, data, weights = NULL, svydesign = NULL, pop_totals = NULL, pop_size = NULL) {

  if (!is.null(svydesign)) {
  model_Frame <- model.frame(formula, data)
  y_nons <- model.response(model_Frame)
  outcome_name <- names(model_Frame)[1]
  mt <- attr(model_Frame, "terms")
  nons_names <- attr(mt, "term.labels") # colnames(get_all_vars(formula, data)) names of variables of nonprobability sample terms(formula, data = data)
  if (all(nons_names %in% colnames(svydesign$variables))) {
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
       outcome_name = outcome_name)

  } else if (!is.null(pop_totals)) { # TODO
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
    if(all(total_names %in% names(pop_totals))) { # pop_totals, pop_means defined such as in `calibrate` function
      pop_totals <- pop_totals[total_names]
    } else {
      warning("Selection and population totals have different names.")
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
  res <- list(
    theta = object$parameters,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )

  attr(res$theta, "glm") <- TRUE
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE
  attr(res, "model")     <- c("glm regression on selection variable")
  res
}

specific_summary_info.nonprobsvy_mi <- function(object,
                                                ...) {

  res <- list(
    beta = object$parameters
  )
  if (object$control$method_outcome == "glm") {
  attr(res$beta, "glm") <- TRUE
  attr(res, "model") <- "glm regression on outcome variable"
  } else if (object$control$method_outcome == "nn") {
    attr(res$beta, "glm") <- FALSE
    attr(res, "model") <- "non-parametric method on outcome variable"
  }
  res
}

specific_summary_info.nonprobsvy_dr <- function(object,
                                                ...) { # TODO for method_outcome equal to nn
  res <- list(
    theta = object$parameters,
    beta  = object$beta,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )

  attr(res$beta,  "glm") <- TRUE
  attr(res$theta, "glm") <- TRUE
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE
  attr(res, "model")     <- c("glm regression on selection variable",
                             "glm regression on outcome variable")
  res
}
