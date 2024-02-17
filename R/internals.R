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
                         ...) { # consider penalty factor here # TO consider for pop_totals/pop_means

  # fit <- glm.fit(x = X, y = y, weights = weights, family = binomial(link = method_selection))
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
  fff <- as.character(formula)
  f <- strsplit(fff[2], "\\s*\\+\\s*")[[1]]
  outcome_formulas <- list()
  if (any(duplicated(f))) {
    warning("No unique names of the outcome variables in formula. The error has been corrected")
    f <- unique(f)
  }
  l <- length(f)
  for (i in 1:l) {
    outcome_formulas[[i]] <- as.formula(paste(f[i], fff[3], sep = " ~ "))
  }

  list(
    f = f,
    outcomes = outcome_formulas,
    l = l
  )
}

