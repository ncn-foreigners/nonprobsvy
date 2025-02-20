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

  model_nons <- nonprob_mi_nn(
    data = X_nons,
    query = X_nons,
    k = control$k,
    treetype = control$treetype,
    searchtype = control$searchtype
  )
  if (is.null(pop_totals)) {
    model_rand <- nonprob_mi_nn(
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
      FUN = function(x) mean(y_nons[x])
      # FUN=function(x) mean(sample_nonprob$short_[x])
    )

    y_nons_pred <- apply(model_nons$nn.idx, 1,
      FUN = function(x) mean(y_nons[x])
      # FUN=function(x) mean(sample_nonprob$short_[x])
    )
  } else {
    model_rand <- nonprob_mi_nn(
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
      FUN = function(x) mean(y_nons[x])
      # FUN=function(x) mean(sample_nonprob$short_[x])
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


nonprob_mi_nn <- function(data,
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


#' "exact" variance estimator for the NN estimator
#' should be renamed with something that refers to variance
#' @noRd
nn_exact <- function(pi_ij,
                     weights_rand,
                     n_nons,
                     y,
                     X_nons,
                     X_rand,
                     k,
                     control,
                     N) {
  loop_size <- 50

  dd <- numeric(length = loop_size)
  for (jj in 1:loop_size) {
    boot_samp <- sample(1:n_nons, size = n_nons, replace = TRUE)
    # boot_samp <- sample(1:n_rand, size = n_rand, replace = TRUE)
    y_nons_b <- y[boot_samp]
    x_nons_b <- X_nons[boot_samp, , drop = FALSE]

    YY <- nonprob_mi_nn(
      data = x_nons_b,
      query = X_rand,
      k = k,
      treetype = control$treetype,
      searchtype = control$searchtype
    )

    dd[jj] <- weighted.mean(
      apply(YY$nn.idx, 1, FUN = function(x) mean(y_nons_b[x])),
      weights_rand
    )
  }
  var(dd)
}
