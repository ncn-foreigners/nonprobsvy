#' @useDynLib nonprobsvy
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats qnorm
#' @importFrom stats as.formula
#' @importFrom stats terms
#' @importFrom stats reformulate
#' @importFrom survey svytotal
#' @importFrom MASS ginv
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @rdname nonprob
#' @noRd
nonprob_ipw <- function(selection,
                        target,
                        data,
                        svydesign,
                        pop_totals,
                        pop_means,
                        pop_size,
                        method_selection,
                        strata,
                        case_weights,
                        na_action,
                        control_selection,
                        control_inference,
                        start_selection,
                        verbose,
                        se,
                        pop_size_fixed,
                        ...) {

  ## passing control parameters
  gee_h_fun <- control_selection$gee_h_fun
  maxit <- control_selection$maxit
  optim_method <- control_selection$optim_method
  var_method <- control_inference$var_method
  est_method <- control_selection$est_method
  num_boot <- control_inference$num_boot
  var_selection <- control_inference$vars_selection
  lambda_control <- control_selection$lambda
  lambda_min <- control_selection$lambda_min
  nlambda <- control_selection$nlambda
  nfolds <- control_selection$nfolds
  eps <- control_selection$epsilon
  rep_type <- control_inference$rep_type

  ## estimation method
  estimation_method <- est_method_ipw(est_method)

  ## multiple dependent variables
  dependents <- paste(selection, collapse = " ")
  outcome <- stats::as.formula(paste(target[2], dependents))
  outcomes <- make_outcomes(outcome)

  model <- make_model_frame(
    formula = outcomes$outcome[[1]],
    data = data,
    svydesign = svydesign,
    pop_totals = pop_totals
  )

  if (!is.null(svydesign)) {

    X_nons <- model$X_nons
    X_rand <- model$X_rand
    nons_names <- model$nons_names
    X <- rbind(X_rand, X_nons)
    R_nons <- rep(1, NROW(X_nons))
    R_rand <- rep(0, NROW(X_rand))
    R <- c(R_rand, R_nons)

    loc_nons <- which(R == 1)
    loc_rand <- which(R == 0)

    n_nons <- NROW(X_nons)
    n_rand <- NROW(X_rand)

    ps_rand <- svydesign$prob
    weights_rand <- 1 / ps_rand

  } else {

    X_nons <- model$X_nons
    nons_names <- model$nons_names
    y_nons <- model$y_nons
    R <- rep(1, NROW(X_nons))
    n_nons <- NROW(X_nons)
    pop_totals <- model$pop_totals

    X_rand <- NULL
    est_ps_rand <- NULL
    est_ps_rand_der <- NULL
    ps_rand <- NULL
    n_rand <- 0
    weights_rand <- NULL
    X <- rbind(model$X_rand, model$X_nons)
  }

  output <- list()
  ys <- list()

  if (!is.null(svydesign)) {
    if (var_selection == TRUE) {
      X_stand <- ncvreg::std(X) # intercept is removed
      prior_weights <- c(weights_rand, case_weights)

      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      inv_link <- method$make_link_inv

      # Cross-validation for variable selection
      cv <- cv_nonprobsvy_rcpp(
        X = X_stand,
        R = R,
        weights_X = prior_weights,
        method_selection = method_selection,
        gee_h_fun = gee_h_fun,
        maxit = maxit,
        eps = eps,
        lambda_min = lambda_min,
        nlambda = nlambda,
        nfolds = control_selection$nfolds,
        penalty = control_selection$penalty,
        a = switch(control_selection$penalty,
                   SCAD = control_selection$a_SCAD,
                   control_selection$a_MCP
        ),
        lambda = lambda_control,
        pop_totals = pop_totals,
        verbose = verbose
      )
      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- cv$theta_selected
      cve_selection <- cv$cv_error

      idx <- c(1, theta_selected + 2) # intercept plus selected variables
      # psel <- length(idx)
      X <- as.matrix(X[, idx])
      # X_design <- cbind(1, Xsel)
      X_nons <- X[loc_nons, , drop = FALSE]
      X_rand <- X[loc_rand, , drop = FALSE]
    }

    selection_model <- estimation_method$estimation_model(
      model = estimation_method$model_selection(
        X = X,
        X_nons = X_nons,
        X_rand = X_rand,
        weights = case_weights,
        weights_rand = weights_rand,
        R = R,
        method_selection = method_selection,
        optim_method = optim_method,
        gee_h_fun = gee_h_fun,
        est_method = est_method,
        maxit = maxit,
        start = start_selection,
        control_selection = control_selection,
        verbose = verbose,
        varcov = TRUE
      ),
      method_selection = method_selection
    )

    theta_hat <- selection_model$theta_hat
    grad <- selection_model$grad
    hess <- selection_model$hess
    var_cov1 <- selection_model$var_cov1
    var_cov2 <- selection_model$var_cov2
    ps_nons <- selection_model$ps_nons
    est_ps_rand <- selection_model$est_ps_rand
    ps_nons_der <- selection_model$ps_nons_der
    est_ps_rand_der <- selection_model$est_ps_rand_der
    theta_standard_errors <- sqrt(diag(selection_model$variance_covariance))

    names(theta_hat) <- names(selection_model$theta_hat) <- colnames(X)
    weights_nons <- 1 / ps_nons
    N <- sum(case_weights * weights_nons)

  } else {
    if (var_selection == TRUE) {
      X_stand <- ncvreg::std(X) # penalizing without an intercept
      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      inv_link <- method$make_link_inv

      # Cross-validation for variable selection
      cv <- cv_nonprobsvy_rcpp(
        X = X_stand,
        R = R,
        weights_X = case_weights,
        method_selection = method_selection,
        gee_h_fun = gee_h_fun,
        maxit = maxit,
        eps = eps,
        lambda_min = lambda_min,
        nlambda = nlambda,
        nfolds = control_selection$nfolds,
        penalty = control_selection$penalty,
        a = switch(control_selection$penalty,
                   SCAD = control_selection$a_SCAD,
                   control_selection$a_MCP
        ),
        lambda = lambda_control,
        pop_totals = pop_totals[-1],
        verbose = verbose
      )

      min <- cv$min
      lambda <- cv$lambda
      theta_selected <- cv$theta_selected
      cve_selection <- cv$cv_error

      idx <- c(1, theta_selected + 2) # intercept + selected variables
      psel <- length(idx)
      X <- as.matrix(X[, idx, drop = FALSE])
      X_nons <- X_nons[, idx, drop = FALSE]
      pop_totals <- model$pop_totals[idx]
    }

    h_object <- theta_h_estimation(
      R = R,
      X = X,
      weights_rand = weights_rand,
      weights = case_weights,
      gee_h_fun = gee_h_fun,
      method_selection = method_selection,
      start = if (is.null(start_selection)) numeric(ncol(X)) else start_selection,
      maxit = maxit,
      nleqslv_method = control_selection$nleqslv_method,
      nleqslv_global = control_selection$nleqslv_global,
      nleqslv_xscalm = control_selection$nleqslv_xscalm,
      pop_totals = pop_totals
    ) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation

    theta_hat <- h_object$theta_h
    hess <- h_object$hess
    grad <- h_object$grad
    names(theta_hat) <- colnames(X)

    method <- switch(method_selection,
                     "logit" = method_ps("logit"),
                     "probit" = method_ps("probit"),
                     "cloglog" = method_ps("cloglog"))

    inv_link <- method$make_link_inv
    dinv_link <- method$make_link_inv_der
    eta_nons <- theta_hat %*% t(X_nons)
    ps_nons <- inv_link(eta_nons)
    ps_nons_der <- dinv_link(eta_nons)
    variance_covariance <- try(solve(-hess), silent = TRUE)

    if (inherits(variance_covariance, "try-error")) {
      if (verbose) message("solve() failed, using ginv() instead.")
      variance_covariance <- MASS::ginv(-hess)
    }

    theta_standard_errors <- sqrt(diag(variance_covariance))
    var_cov1 <- method$variance_covariance1
    var_cov2 <- method$variance_covariance2
    df_residual <- NROW(X_nons) - length(theta_hat)
    weights_nons <- 1 / ps_nons
    N <- sum(case_weights * weights_nons)
    residuals <- as.vector(rep(1, n_nons) - ps_nons)

    selection_model <- list(
      theta_hat = theta_hat,
      hess = hess,
      grad = grad,
      ps_nons = ps_nons,
      est_ps_rand = est_ps_rand,
      ps_nons_der = ps_nons_der,
      est_ps_rand_der = est_ps_rand_der,
      variance_covariance = variance_covariance,
      var_cov1 = var_cov1,
      var_cov2 = var_cov2,
      df_residual = df_residual,
      residuals = residuals,
      log_likelihood = NA
    )
    #######################################
  }

  mu_hats <- numeric(length = outcomes$l)
  for (o in 1:outcomes$l) {
    if (is.null(pop_totals)) {
      y_nons <- make_model_frame(
        formula = outcomes$outcome[[o]],
        data = data,
        svydesign = svydesign,
        weights = case_weights
      )$y_nons
    } else {
      y_nons <- make_model_frame(
        formula = outcomes$outcome[[o]],
        data = data,
        pop_totals = pop_totals,
        weights = case_weights
      )$y_nons
    }
    ys[[o]] <- as.numeric(y_nons)
    mu_hats[o] <- mu_hatIPW(
      y = y_nons,
      weights = case_weights,
      weights_nons = weights_nons,
      N = ifelse(is.null(pop_size), N, pop_size)
    )
  }

  if (se) {
    confidence_interval <- list()
    SE_values <- list()
    if (var_method == "analytic") {
      var_nonprob <- numeric(length = outcomes$l)
      var_prob <- numeric(length = outcomes$l)
      var <- numeric(length = outcomes$l)
      se_nonprob <- numeric(length = outcomes$l)
      se_prob <- numeric(length = outcomes$l)
      for (k in 1:outcomes$l) {
        var_obj <- internal_varIPW(
          svydesign = svydesign,
          X_nons = X_nons,
          X_rand = X_rand,
          y_nons = ys[[o]],
          weights = case_weights,
          ps_nons = ps_nons,
          mu_hat = mu_hats[o],
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
          gee_h_fun = gee_h_fun,
          var_cov1 = var_cov1,
          var_cov2 = var_cov2,
          verbose = verbose
        )

        var_nonprob[k] <- var_obj$var_nonprob
        var_prob[k] <- var_obj$var_prob
        var[k] <- var_obj$var
        se_nonprob[k] <- sqrt(var_nonprob[o])
        se_prob[k] <- sqrt(var_prob[o])
        SE_values[[k]] <- data.frame(prob = se_prob[o], nonprob = se_nonprob[o])
      }
    } else if (var_method == "bootstrap") { # TODO add ys, mu_hats instead of y_nons,
      if (control_inference$cores > 1) {
        boot_obj <- boot_ipw_multicore(
          X_rand = X_rand,
          X_nons = X_nons,
          svydesign = svydesign,
          ys = ys, #
          num_boot = num_boot,
          case_weights = case_weights,
          weights_rand = weights_rand,
          R = R,
          theta_hat = theta_hat,
          mu_hats = mu_hats, #
          method_selection = method_selection,
          start_selection = start_selection,
          n_nons = n_nons,
          n_rand = n_rand,
          optim_method = optim_method,
          est_method = est_method,
          gee_h_fun = gee_h_fun,
          maxit = maxit,
          pop_size = pop_size,
          pop_totals = pop_totals,
          control_selection = control_selection,
          control_inference = control_inference,
          cores = control_inference$cores,
          verbose = verbose
        )
      } else {
        boot_obj <- boot_ipw(
          X_rand = X_rand,
          X_nons = X_nons,
          svydesign = svydesign,
          ys = ys, #
          num_boot = num_boot,
          case_weights = case_weights,
          weights_rand = weights_rand,
          R = R,
          theta_hat = theta_hat,
          mu_hats = mu_hats, #
          method_selection = method_selection,
          start_selection = start_selection,
          n_nons = n_nons,
          n_rand = n_rand,
          optim_method = optim_method,
          est_method = est_method,
          gee_h_fun = gee_h_fun,
          maxit = maxit,
          pop_size = pop_size,
          pop_totals = pop_totals,
          control_selection = control_selection,
          control_inference = control_inference,
          verbose = verbose
        )
      }
      var <- boot_obj$var
      # mu_hat <- boot_obj$mu
      for (o in 1:outcomes$l) {
        SE_values[[o]] <- data.frame(nonprob = NA, prob = NA)
      }
    } else {
      stop("Invalid `var_method` for the variance estimation.")
    }
    SE <- sqrt(var)
    z <- stats::qnorm(1 - control_inference$alpha / 2)
    # confidence interval based on the normal approximation
    for (o in 1:outcomes$l) {
      confidence_interval[[o]] <- data.frame(lower_bound = mu_hats[o] - z * SE[o],
                                             upper_bound = mu_hats[o] + z * SE[o])
    }
  } else {
    confidence_interval <- NULL
    SE_values <- NULL

    for (o in 1:outcomes$l) {
      SE <- NA
      confidence_interval[[o]] <- data.frame(lower_bound = NA, upper_bound = NA)
      SE_values[[o]] <- data.frame(nonprob = NA, prob = NA)
    }
  }
  for (o in 1:outcomes$l) {
    output[[o]] <- data.frame(mean = mu_hats[o], SE = SE[o])
  }

  parameters <- matrix(c(theta_hat, theta_standard_errors),
    ncol = 2,
    dimnames = list(
      names(theta_hat),
      c("Estimate", "Std. Error")
    )
  )
  weights_summary <- summary(as.vector(weights_nons))
  prop_scores <- c(ps_nons, est_ps_rand)
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f
  if (is.null(pop_size)) pop_size <- N # estimated pop_size
  names(pop_size) <- "pop_size"
  names(ys) <- all.vars(outcome[[2]])


  boot_sample <- if (se == T & control_inference$var_method == "bootstrap" & control_inference$keep_boot) {
    boot_obj$stat
  } else {
    NULL
  }
  if (!is.null(boot_sample) & is.matrix(boot_sample)) colnames(boot_sample) <- names(ys)


  selection_list <- list(
    selection_model = selection_model,
    coefficients = selection_model$theta_hat,
    std_err = theta_standard_errors,
    residuals = selection_model$residuals,
    variance = selection_model$variance,
    fitted_values = prop_scores,
    link = selection_model$method,
    linear_predictors = selection_model$eta,
    aic = selection_model$aic,
    ipw_weights = as.vector(weights_nons),
    case_weights = case_weights,
    pop_totals = pop_totals,
    formula = selection,
    df_residual = selection_model$df_residual,
    log_likelihood = selection_model$log_likelihood,
    method_selection = method_selection,
    hessian = hess,
    gradient = grad,
    method = est_method,
    prob_der = ps_nons_der,
    prob_rand = ps_rand,
    prob_rand_est = est_ps_rand,
    prob_rand_est_der = est_ps_rand_der,
    gee_h_fun = gee_h_fun,
    cve = if (control_inference$vars_selection == TRUE) {
      cve_selection
    } else {
      NULL
    }
  )

  structure(
    list(
      data = data,
      X = rbind(X_rand, X_nons),
      y = ys,
      R = R,
      ps_scores = prop_scores,
      case_weights = case_weights,
      ipw_weights = as.vector(weights_nons),
      control = list(
        control_selection = control_selection,
        control_outcome = NULL,
        control_inference = control_inference
      ),
      output = output,
      SE = SE_values,
      confidence_interval = confidence_interval,
      nonprob_size = n_nons,
      prob_size = n_rand,
      pop_size = pop_size,
      pop_size_fixed = pop_size_fixed,
      pop_totals = pop_totals,
      pop_means = pop_means,
      outcome = NULL,
      selection = selection_list,
      boot_sample = boot_sample,
      svydesign = if (is.null(pop_totals)) svydesign else NULL,
      ys_rand_pred = NULL,
      ys_nons_pred = NULL,
      ys_resid = NULL
    ),
    class = "nonprob"
  )
}
