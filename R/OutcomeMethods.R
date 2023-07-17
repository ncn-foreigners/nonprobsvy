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
                model_frame) {

  # Estimation for outcome model
  model_out <- internal_outcome(outcome = outcome,
                                data = data,
                                weights = weights,
                                family_outcome = family_outcome)

  model_nons_coefs <- model_out$glm$coefficients
  parameters <- model_out$glm_summary$coefficients

  # y_rand_pred <- as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample # consider predict function
  # print(names(attr(X_rand, "contrasts")))
  # print(attr(X_rand, "contrasts"))
  # print(colnames(X_rand))
  # print(dim(model_frame))
  y_rand_pred <- predict.glm(model_out$glm, newdata = model_frame, type = "response")
  # if(family_outcome == "binomial") {
  #   y_rand_pred <- ifelse(y_rand_pred >= .5, 1, 0)
  # }
  y_nons_pred <- model_out$glm$fitted.values #as.numeric(X_nons %*% model_nons_coefs)
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

  # idx <- model_rand$nn.idx
  # y_rand_pred <- y_nons[idx] TODO without loop

  for (i in 1:n_rand) {
    idx <- model_rand$nn.idx[i,]
    y_rand_pred[i] <- weighted.mean(y_nons[idx], weights[idx])
  }

  for (i in 1:n_nons) {
    idx <- model_nons$nn.idx[i,]
    y_nons_pred[i] <- weighted.mean(y_nons[idx], weights[idx])
  }
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
                           #weights = weights,
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

