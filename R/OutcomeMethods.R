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

  if(is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }
    if (vars_selection == FALSE) {
      # Estimation for outcome model
      model_out <- internal_outcome(outcome = outcome,
                                    data = data,
                                    weights = weights,
                                    family_outcome = family_outcome,
                                    start_outcome = start_outcome)

      model_nons_coefs <- model_out$glm$coefficients
      parameters <- model_out$glm_summary$coefficients

      if (is.null(pop_totals)) {
        y_rand_pred <- stats::predict.glm(model_out$glm, newdata = model_frame, type = "response")
      } else {
        eta <- pop_totals %*% model_nons_coefs / pop_totals[1]
        y_rand_pred <- family_nonprobsvy$mu(eta)
      }
      y_nons_pred <- model_out$glm$fitted.values
    } else {
      model <- stats::glm.fit(x = X_nons,
                              y = y_nons,
                              weights = weights,
                              family = get_method(family_outcome),
                              start = start_outcome,
                              control = list(control$epsilon,
                                              control$maxit,
                                              control$trace),
                              intercept = FALSE)
      model_summ <- stats::summary.glm(model)
      parameters <- model_summ$coefficients
      model_nons_coefs <- model$coefficients
      if (is.null(pop_totals)) {
        eta <- X_rand %*% model_nons_coefs
      } else {
        eta <- pop_totals %*% model_nons_coefs / pop_totals[1]
      }
      y_rand_pred <- family_nonprobsvy$mu(eta)
      y_nons_pred <- model$fitted.values

      model_out <- list(glm = model,
                        glm_summary = model_summ)
    }
  model_out$glm$std_err <- parameters[,2]
  names(model_out$glm$std_err) <- names(model_out$glm$coefficients)

  list(model = model_out$glm,
       y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       parameters = parameters)
}

nn_nonprobsvy <- function(outcome,
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
                          pop_totals,
                          model_frame = NULL,
                          start_outcome = NULL) { # TODO consider add data standardization before modelling
  model_nons <- nonprobMI_nn(data = X_nons,
                             query = X_nons,
                             k = control$k,
                             treetype = control$treetype,
                             searchtype = control$searchtype)
  if (is.null(pop_totals)) {
    model_rand <- nonprobMI_nn(data = X_nons,
                               query = X_rand,
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
  } else {
    model_rand <- nonprobMI_nn(data = X_nons,
                               query = t(as.matrix(pop_totals / pop_totals[1])),
                               k = control$k,
                               treetype = control$treetype,
                               searchtype = control$searchtype)
    y_rand_pred <- vector(mode = "numeric", length = 1)
    y_nons_pred <- vector(mode = "numeric", length = n_nons)
    parameters <- "Non-parametric method for outcome model"

    y_rand_pred <- mean(y_nons[model_rand$nn.idx])
    y_nons_pred <- apply(model_nons$nn.idx, 1,
                         FUN=\(x) mean(y_nons[x])
                         #FUN=\(x) mean(sample_nonprob$short_[x])
    )
  }

  model_out <- list(model_nons = model_nons,
                    model_rand = model_rand)
  list(model = model_out,
       y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       parameters = parameters)
}

pmm_nonprobsvy <- function(outcome,
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
                           vars_selection,
                           pop_totals,
                           model_frame) {

  glm_object <- glm_nonprobsvy(outcome,
                               data,
                               weights,
                               family_outcome,
                               start_outcome = start_outcome,
                               X_nons,
                               y_nons,
                               X_rand,
                               control,
                               n_nons,
                               n_rand,
                               model_frame,
                               vars_selection,
                               pop_totals)

  # This is commented now because it is not needed
  # model_nons <- nonprobMI_nn(data = glm_object$y_nons_pred,
  #                            query = glm_object$y_nons_pred,
  #                            k = control$k,
  #                            treetype = control$treetype,
  #                            searchtype = control$searchtype)
  #
  # y_nons_pred <- apply(model_nons$nn.idx, 1,
  #                      FUN=\(x) mean(y_nons[x])
  #                      #FUN=\(x) mean(sample_nonprob$short_[x])
  # )

  switch (control$predictive_match,
    { # 1
      if (is.null(pop_totals)) {
        model_rand <- nonprobMI_nn(data = y_nons,
                                   query = glm_object$y_rand_pred,
                                   k = control$k,
                                   treetype = control$treetype,
                                   searchtype = control$searchtype)

        y_rand_pred <- apply(model_rand$nn.idx, 1,
                             FUN=\(x) mean(y_nons[x])
                             #FUN=\(x) mean(sample_nonprob$short_[x])
        )
      } else {
        model_rand <- nonprobMI_nn(data = y_nons,
                                   query = glm_object$y_rand_pred,
                                   k = control$k,
                                   treetype = control$treetype,
                                   searchtype = control$searchtype)
        y_rand_pred <- mean(y_nons[model_rand$nn.idx])
      }
    },
    { # 2
      if (is.null(pop_totals)) {
        model_rand <- nonprobMI_nn(data = glm_object$y_nons_pred,
                                   query = glm_object$y_rand_pred,
                                   k = control$k,
                                   treetype = control$treetype,
                                   searchtype = control$searchtype)

        y_rand_pred <- apply(model_rand$nn.idx, 1,
                             FUN=\(x) mean(y_nons[x])
                             #FUN=\(x) mean(sample_nonprob$short_[x])
        )
      } else {
        model_rand <- nonprobMI_nn(data = glm_object$y_nons_pred,
                                   query = glm_object$y_rand_pred,
                                   k = control$k,
                                   treetype = control$treetype,
                                   searchtype = control$searchtype)
        y_rand_pred <- mean(y_nons[model_rand$nn.idx])
      }
    }
  )

  model_out <- list(#model_nons = model_nons,
                    model_rand = model_rand)
  list(model = model_out,
       y_rand_pred = y_rand_pred,
       #y_nons_pred = y_nons_pred,
       parameters = glm_object$parameters)
}

nonprobMI_fit <- function(outcome,
                          data,
                          weights,
                          svydesign,
                          family_outcome,
                          start,
                          control_outcome = controlOut(),
                          verbose,
                          model,
                          x,
                          y) {

  family <- family_outcome

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  data$weights <- weights # TODO just for now, find more efficient way
  model_nons <- stats::glm(formula = outcome,
                           data = data,
                           weights = weights,
                           family = family,
                           start = start,
                           control = list(control_outcome$epsilon,
                                          control_outcome$maxit,
                                          control_outcome$trace))

  model_nons
}

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

