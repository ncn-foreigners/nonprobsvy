#' @importFrom stats predict.glm
#' @importFrom stats glm.fit
#' @importFrom stats summary.glm
glm <- function(outcome,
                data,
                weights,
                family_outcome,
                X_nons,
                y_nons,
                X_rand,
                control,
                n_nons,
                n_rand,
                model_frame,
                vars_selection) {

  if (vars_selection == FALSE) {
    # Estimation for outcome model
    model_out <- internal_outcome(outcome = outcome,
                                  data = data,
                                  weights = weights,
                                  family_outcome = family_outcome)

    model_nons_coefs <- model_out$glm$coefficients
    parameters <- model_out$glm_summary$coefficients

    y_rand_pred <- stats::predict.glm(model_out$glm, newdata = model_frame, type = "response")
    y_nons_pred <- model_out$glm$fitted.values
  } else {
    if(is.character(family_outcome)) {
      family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
      family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
      family_nonprobsvy <- family_nonprobsvy()
    }
    model <- stats::glm.fit(x = X_nons,
                             y = y_nons,
                             weights = weights,
                             family = get_method(family_outcome),
                             control = list(control$epsilon,
                                            control$maxit,
                                            control$trace),
                             intercept = FALSE)
    model_summ <- stats::summary.glm(model)
    parameters <- model_summ$coefficients
    model_nons_coefs <- model$coefficients
    eta <- X_rand %*% model_nons_coefs

    y_rand_pred <- family_nonprobsvy$mu(eta)
    y_nons_pred <- model$fitted.values

    model_out <- list(glm = model,
                      glm_summary = model_summ)
  }
  list(model = model_out,
       y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       parameters = parameters)
}


nn <- function(outcome,
               data,
               weights,
               family_outcome,
               X_nons,
               y_nons,
               X_rand,
               control,
               n_nons,
               n_rand,
               vars_selection,
               model_frame = NULL) {

  model_rand <- nonprobMI_nn(data = X_nons,
                             query = X_rand,
                             k = control$k,
                             treetype = control$treetype,
                             searchtype = control$searchtype)
  model_nons <- nonprobMI_nn(data = X_nons,
                             query = X_nons,
                             k = control$k,
                             treetype = control$treetype,
                             searchtype = control$searchtype)
  y_rand_pred <- vector(mode = "numeric", length = n_rand)
  y_nons_pred <- vector(mode = "numeric", length = n_nons)
  parameters <- "Non-parametric method for outcome model"

  y_rand_pred <- apply(model_rand$nn.idx, 1,
                      FUN=\(x) mean(y_nons[x])
                      #FUN=\(x) mean(sample_nonprob$short_[x])
  )

  y_nons_pred <- apply(model_nons$nn.idx, 1,
                      FUN=\(x) mean(y_nons[x])
                      #FUN=\(x) mean(sample_nonprob$short_[x])
  )

  model_out <- list(model_nons = model_nons,
                    model_rand = model_rand)
  list(model = model_out,
       y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       parameters = parameters)
}

#' nonprobMI_fit
#
#' nonprobMI_fit: Function for outcome variable estimation based on nonprobability sample and using model based approach
#'
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param control_outcome - a
#' @param start - a
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param verbose - a
#' @param model - a
#' @param x - a
#' @param y - a
#'


nonprobMI_fit <- function(outcome,
                          data,
                          weights,
                          svydesign,
                          family_outcome,
                          control_outcome = controlOut(),
                          start = NULL,
                          verbose,
                          model,
                          x,
                          y) { # TODO problem with weights argument

  family <- family_outcome

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  model_nons <- stats::glm(formula = outcome,
                           family = family,
                           data = data,
                           # weights = weights, invalid type (closure) for variable '(weights)' TOFIX
                           start = start,
                           control = list(control_outcome$epsilon,
                                          control_outcome$maxit,
                                          control_outcome$trace))

  model_nons
}



#' nonprobMI_nn
#
#' nonprobMI_nn: Function for outcome variable estimation based on nonprobability sample and using predictive mean matching
#'
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param query - a
#' @param k - a
#' @param treetype - a
#' @param searchtype - a
#' @param radius - a
#' @param eps - a

nonprobMI_nn <- function(data,
                         query,
                         k,
                         treetype,
                         searchtype,
                         radius = 0,
                         eps = 0) {


  model_nn <- RANN::nn2(data = data,
                        query = query,
                        k = k,
                        treetype = treetype,
                        searchtype = searchtype,
                        radius = radius,
                        eps = eps)
  model_nn
}

