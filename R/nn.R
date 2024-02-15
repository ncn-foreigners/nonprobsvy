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

  model_nons <- nonprobMI_nn(
    data = X_nons,
    query = X_nons,
    k = control$k,
    treetype = control$treetype,
    searchtype = control$searchtype
  )
  if (is.null(pop_totals)) {
    model_rand <- nonprobMI_nn(
      data = X_nons,
      query = X_rand,
      k = control$k,
      treetype = control$treetype,
      searchtype = control$searchtype
    )
    y_rand_pred <- vector(mode = "numeric", length = n_rand)
    y_nons_pred <- vector(mode = "numeric", length = n_nons)
    parameters <- "Non-parametric method for outcome model"

    y_rand_pred <- apply(model_rand$nn.idx, 1,
      FUN = \(x) mean(y_nons[x])
      # FUN=\(x) mean(sample_nonprob$short_[x])
    )

    y_nons_pred <- apply(model_nons$nn.idx, 1,
      FUN = \(x) mean(y_nons[x])
      # FUN=\(x) mean(sample_nonprob$short_[x])
    )
  } else {
    model_rand <- nonprobMI_nn(
      data = X_nons,
      query = t(as.matrix(pop_totals / pop_totals[1])),
      k = control$k,
      treetype = control$treetype,
      searchtype = control$searchtype
    )
    y_rand_pred <- vector(mode = "numeric", length = 1)
    y_nons_pred <- vector(mode = "numeric", length = n_nons)
    parameters <- "Non-parametric method for outcome model"

    y_rand_pred <- mean(y_nons[model_rand$nn.idx])
    y_nons_pred <- apply(model_nons$nn.idx, 1,
      FUN = \(x) mean(y_nons[x])
      # FUN=\(x) mean(sample_nonprob$short_[x])
    )
  }

  model_out <- list(
    model_nons = model_nons,
    model_rand = model_rand
  )
  list(
    model = model_out,
    y_rand_pred = y_rand_pred,
    y_nons_pred = y_nons_pred,
    parameters = parameters
  )
}


nonprobMI_nn <- function(data,
                         query,
                         k,
                         treetype,
                         searchtype,
                         radius = 0,
                         eps = 0) {
  model_nn <- RANN::nn2(
    data = data,
    query = query,
    k = k,
    treetype = treetype,
    searchtype = searchtype,
    radius = radius,
    eps = eps
  )
  model_nn
}
