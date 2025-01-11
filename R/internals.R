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
#' @importFrom stats cov
#' @importFrom stats var
#' @importFrom stats predict

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
                               verbose,
                               bias_correction = FALSE,
                               varcov = FALSE,
                               ...) {
  if (bias_correction == TRUE) est_method <- "mm"
  estimation_method <- get_method(est_method)
  estimation_method$model_selection(
    X = X,
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
    verbose = verbose,
    ...
  )
}

# Outcome model object
internal_outcome <- function(outcome,
                             data,
                             weights,
                             family_outcome,
                             start_outcome) {
  # estimation
  model_nons <- nonprobMI_fit(
    outcome = outcome,
    data = data,
    weights = weights,
    family_outcome = family_outcome,
    start = start_outcome
  )
  model_nons_summary <- summary(model_nons)

  list(
    glm = model_nons,
    glm_summary = model_nons_summary
  )
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
                         ...) {
  # TODO
  # consider penalty factor here
  # consider for pop_totals/pop_means
  if (is.null(pop_totals)) {
    fit <- stats::glm(y ~ 1,
      weights = weights,
      family = binomial(link = method_selection)
    )

    n <- length(y)
    p <- ncol(X)
    w <- fit$weights
    r <- as.matrix(stats::residuals(fit, "working") * w)
    zmax <- max(crossprod(X, r)) / n
    lambda_max <- zmax / alpha
  } else {
    lambda_max <- .1
  }
  if (log_lambda) { # lambda sequence on log-scale
    if (lambda_min == 0) {
      lambda <- c(exp(seq(log(lambda_max), log(.001 * lambda_max), length = nlambda - 1)), 0)
    } else {
      lambda <- exp(seq(log(lambda_max), log(lambda_min * lambda_max), length = nlambda))
    }
  } else { # lambda sequence on linear-scale
    if (lambda_min == 0) {
      lambda <- c(seq(lambda_max, 0.001 * lambda_max, length = nlambda - 1), 0)
    } else {
      lambda <- seq(lambda_max, lambda_min * lambda_max, length = nlambda)
    }
  }
  lambda
}

start_fit <- function(X,
                      R,
                      weights,
                      weights_rand,
                      method_selection,
                      control_selection = controlSel()) {
  weights_to_glm <- c(weights_rand, weights)
  start_model <- stats::glm.fit(
    x = X, # glm model for initial values in propensity score estimation
    y = R,
    weights = weights_to_glm, # to fix
    family = binomial(link = method_selection),
    control = list(
      epsilon = control_selection$epsilon,
      maxit = control_selection$maxit,
      trace = control_selection$trace
    )
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

ff <- function(formula) {
  formula_string <- paste(deparse(formula), collapse = " ")
  formula_parts <- strsplit(formula_string, "~")[[1]]
  if (length(formula_parts) != 2) {
    stop("The formula must contain exactly one '~' operator.")
  }

  lhs <- trimws(formula_parts[1])
  rhs <- trimws(formula_parts[2])

  dependent_vars <- strsplit(lhs, "\\s*\\+\\s*")[[1]]
  independent_vars <- strsplit(rhs, "\\s*\\+\\s*")[[1]]

  if (any(duplicated(dependent_vars))) {
    warning("Duplicate dependent variable names detected. They have been made unique.")
    dependent_vars <- unique(dependent_vars)
  }
  outcome_formulas <- vector("list", length(dependent_vars))
  l <- length(dependent_vars)
  for (i in seq_along(dependent_vars)) {
    outcome_formulas[[i]] <- reformulate(termlabels = independent_vars, response = dependent_vars[i])
  }
  list(
    f = dependent_vars,
    outcomes = outcome_formulas,
    l = l
  )
}

mu_hatDR <- function(y,
                     y_nons,
                     y_rand,
                     weights,
                     weights_nons,
                     weights_rand,
                     N_nons,
                     N_rand) {
  correction_term <- sum(weights * weights_nons * (y - y_nons)) / N_nons
  probability_estimate <- sum(weights_rand * y_rand) / N_rand
  correction_term + probability_estimate
}

mu_hatIPW <- function(y,
                      weights,
                      weights_nons,
                      N) {
  mu_hat <- sum(weights * weights_nons * y) / N
  mu_hat
}

nonprobMI_fit <- function(outcome,
                          data,
                          weights,
                          svydesign = NULL,
                          family_outcome = "gaussian",
                          start = NULL,
                          control_outcome = controlOut(),
                          verbose = FALSE,
                          model = TRUE,
                          x = FALSE,
                          y = FALSE) {
  # Process family specification
  family <- process_family(family_outcome)

  # Process control parameters
  control_list <- list(
    epsilon = control_outcome$epsilon,
    maxit = control_outcome$maxit,
    trace = control_outcome$trace
  )

  # Create model environment to avoid modifying original data
  model_data <- data
  model_data$weights <- weights

  # Fit the model
  tryCatch({
    model_fit <- stats::glm(
      formula = outcome,
      data = model_data,
      weights = weights,
      family = family,
      start = start,
      control = control_list,
      model = model,
      x = x,
      y = y
    )

    if (verbose) {
      cat("Model fitting completed:\n")
      cat("Convergence status:", ifelse(model_fit$converged, "converged", "not converged"), "\n")
      cat("Number of iterations:", model_fit$iter, "\n")
    }

    return(model_fit)

  }, error = function(e) {
    stop("Error in model fitting: ", e$message)
  })
}

process_family <- function(family_spec) {
  if (is.character(family_spec)) {
    family <- get(family_spec, mode = "function", envir = parent.frame())
  } else if (is.function(family_spec)) {
    family <- family_spec()
  } else if (inherits(family_spec, "family")) {
    family <- family_spec
  } else {
    stop("Invalid family specification")
  }
  return(family)
}
