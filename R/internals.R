##' @title Setup lambda parameter for the {ncvreg} package
##' @param X a `matrix` of auxiliary variables
##' @param y a `vector` of target variables
##'
##' @noRd
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



##' @title Function that create formulas for multiple outcomes
##'
##' @param formula a formula specifying the outcome ~ auxiliary variables (e.g. `y1 + y2 ~ x1 + x2`).
##'
##' @description
##' The function takes an outcome formula such as `y1 + y2 ~ x1 + x2` and creates
##' separate models i.e. `y1 ~ x1 + x2` and `y2 ~ x1 + x2`
##' @noRd
make_outcomes <- function(formula) {
  formula_string <- paste(deparse(formula), collapse = " ")
  formula_parts <- strsplit(formula_string, "~")[[1]]
  if (length(formula_parts) != 2) {
    stop("The `formula` must contain exactly one '~' operator.")
  }

  lhs <- trimws(formula_parts[1])
  rhs <- trimws(formula_parts[2])

  dependent_vars <- strsplit(lhs, "\\s*\\+\\s*")[[1]]
  independent_vars <- strsplit(rhs, "\\s*\\+\\s*")[[1]]

  if (any(duplicated(dependent_vars))) {
    warning("Duplicate dependent variable names have been detected. They have been made unique.")
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

mu_hatDR <- function(y_hat,
                     y_resid,
                     weights,
                     weights_nons,
                     N_nons) {

  colSums(weights * weights_nons * y_resid) / N_nons + y_hat
}

## replace with weighted.mean -- maybe it can be removed
mu_hatIPW <- function(y,
                      weights,
                      weights_nons,
                      N) {
  mu_hat <- sum(weights * weights_nons * y) / N
  mu_hat
}



merge_formulas <- function(outcome, selection) {

  lhs1 <- formula.tools::lhs(outcome)  # Left hand side of formula1
  rhs1 <- formula.tools::rhs(outcome)  # Right hand side of formula1
  rhs2 <- formula.tools::rhs(selection)  # Right hand side of formula2

  # Combine right hand sides (unique terms)
  combined_rhs <- unique(c(all.vars(rhs1), all.vars(rhs2)))
  rhs_terms <- paste(combined_rhs, collapse = " + ")

  outcome <- as.formula(paste(deparse(lhs1), "~", rhs_terms))
  selection <- as.formula(paste("~", rhs_terms))

  return(list(outcome, selection))
}
