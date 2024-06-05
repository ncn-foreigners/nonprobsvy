# Internal functions for mass imputation models
#' @importFrom stats predict.glm
#' @importFrom stats glm.fit
#' @importFrom stats summary.glm
glm_nonprobsvy <- function(outcome,
                           data,
                           weights,
                           family_outcome,
                           start_outcome,
                           X_nons,
                           y_nons,
                           X_rand,
                           control,
                           n_nons,
                           n_rand,
                           model_frame,
                           vars_selection,
                           pop_totals) {
  if (is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
  if (vars_selection == FALSE) {
    # Estimation for outcome model
    model_out <- internal_outcome(
      outcome = outcome,
      data = data,
      weights = weights,
      family_outcome = family_outcome,
      start_outcome = start_outcome
    )

    model_nons_coefs <- model_out$glm$coefficients
    parameters <- model_out$glm_summary$coefficients

    if (is.null(pop_totals)) {
      # print(head(model_frame))
      # stop("123")
      y_rand_pred <- stats::predict.glm(model_out$glm, newdata = model_frame, type = "response")
    } else {
      eta <- pop_totals %*% model_nons_coefs / pop_totals[1]
      y_rand_pred <- family_nonprobsvy$linkinv(eta)
    }
    y_nons_pred <- model_out$glm$fitted.values
  } else {
    model <- stats::glm.fit(
      x = X_nons,
      y = y_nons,
      weights = weights,
      family = get_method(family_outcome),
      start = start_outcome,
      control = list(
        control$epsilon,
        control$maxit,
        control$trace
      ),
      intercept = FALSE
    )
    model_summ <- stats::summary.glm(model)
    parameters <- model_summ$coefficients
    model_nons_coefs <- model$coefficients
    if (is.null(pop_totals)) {
      eta <- X_rand %*% model_nons_coefs
    } else {
      eta <- pop_totals %*% model_nons_coefs / pop_totals[1]
    }
    y_rand_pred <- family_nonprobsvy$linkinv(eta)
    y_nons_pred <- model$fitted.values

    # artificial formula and data created for pmm_extact_se
    predictors <- colnames(X_nons)[-1]
    outcome_name <- names(model_frame)[1]
    formula_str <- paste(outcome_name, "~", paste(predictors, collapse = " + "))
    formula <- as.formula(formula_str)
    model$formula <- formula
    model$data <- as.data.frame(cbind(y_nons, X_nons[,-1, drop = FALSE]))
    colnames(model$data) <- c(outcome_name, predictors)

    model_out <- list(
      glm = model,
      glm_summary = model_summ
    )
  }
  model_out$glm$std_err <- parameters[, 2]
  names(model_out$glm$std_err) <- names(model_out$glm$coefficients)

  list(
    model = model_out$glm,
    y_rand_pred = y_rand_pred,
    y_nons_pred = y_nons_pred,
    parameters = parameters
  )
}
