#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom stats loess.control
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
  glm_object <- switch(control$pmm_reg_engine,
    "glm" = glm_nonprobsvy(
      outcome,
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
      pop_totals
    ),
    "loess" = {
      # doesn't accept weights
      mm <- stats::loess(
        outcome,
        data,
        span = .2,
        control = stats::loess.control(surface = "direct")
      )
      mm$data <- data
      mm$formula <- outcome

      list(
        model = mm,
        y_rand_pred = predict(mm, newdata = model_frame),
        y_nons_pred = predict(mm),
        parameters = NULL
      )
    }
  )

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

  # add protection for very low values in weighting
  switch(control$predictive_match,
    { # 1
      if (is.null(pop_totals)) {
        model_rand <- nonprobMI_nn(
          data = y_nons,
          query = glm_object$y_rand_pred,
          k = control$k,
          treetype = control$treetype,
          searchtype = control$searchtype
        )

        switch(control$pmm_weights,
          "none" = {
            y_rand_pred <- apply(model_rand$nn.idx, 1,
              FUN = \(x) mean(y_nons[x])
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )
          },
          "prop_dist" = {
            # TODO:: these weights will need to be saved for variance estimation
            y_rand_pred <- sapply(1:NROW(model_rand$nn.idx),
              FUN = \(x) weighted.mean(y_nons[model_rand$nn.idx[x, ]],
                w = 1 / model_rand$nn.dist[x, ]
              )
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )
          }
        )
      } else {
        # I'm not touching this
        model_rand <- nonprobMI_nn(
          data = y_nons,
          query = glm_object$y_rand_pred,
          k = control$k,
          treetype = control$treetype,
          searchtype = control$searchtype
        )
        y_rand_pred <- mean(y_nons[model_rand$nn.idx])
      }
    },
    { # 2
      if (is.null(pop_totals)) {
        model_rand <- nonprobMI_nn(
          data = glm_object$y_nons_pred,
          query = glm_object$y_rand_pred,
          k = control$k,
          treetype = control$treetype,
          searchtype = control$searchtype
        )

        y_rand_pred <- apply(model_rand$nn.idx, 1,
          FUN = \(x) mean(y_nons[x])
          # FUN=\(x) mean(sample_nonprob$short_[x])
        )

        switch(control$pmm_weights,
          "none" = {
            y_rand_pred <- apply(model_rand$nn.idx, 1,
              FUN = \(x) mean(y_nons[x])
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )
          },
          "prop_dist" = {
            # TODO:: these weights will need to be saved for variance estimation
            y_rand_pred <- sapply(1:NROW(model_rand$nn.idx),
              FUN = \(x) weighted.mean(y_nons[model_rand$nn.idx[x, ]],
                w = 1 / model_rand$nn.dist[x, ]
              )
              # FUN=\(x) mean(sample_nonprob$short_[x])
            )
          }
        )
      } else {
        # I'm not touching this
        model_rand <- nonprobMI_nn(
          data = glm_object$y_nons_pred,
          query = glm_object$y_rand_pred,
          k = control$k,
          treetype = control$treetype,
          searchtype = control$searchtype
        )
        y_rand_pred <- mean(y_nons[model_rand$nn.idx])
      }
    }
  )

  model_out <- list(
    # model_nons = model_nons,
    model_rand = model_rand,
    glm_object = glm_object$model
  )
  attr(model_out, "method") <- "pmm"
  list(
    model = model_out,
    y_rand_pred = y_rand_pred,
    # y_nons_pred = y_nons_pred,
    parameters = glm_object$parameters
  )
}



pmm_exact <- function(pi_ij,
                      weights_rand,
                      n_nons,
                      y,
                      pmm_reg_engine,
                      stats,
                      glm,
                      model_obj,
                      svydesign,
                      predictive_match,
                      k,
                      N) {
  # if (isTRUE("ppsmat" %in% class(pi_ij))) {
  #   pi_ij <- pi_ij$pij
  # }
  # # if (!is.null(svydesign$dcheck[[1]]$dcheck)) {
  # #   pi_ij <- svydesign$dcheck[[1]]$dcheck
  # # }
  # if (is.null(pi_ij)) {
  #   pi_ij <- outer(1 / weights_rand, 1 / weights_rand) * (
  #     1 - outer(1 - 1 / weights_rand, 1 - 1 / weights_rand) /
  #       sum(1 - 1 / weights_rand))
  # }
  # # if (!is.matrix(pi_ij)) {
  # #
  # # }
  # add variable for loop size to control
  loop_size <- 50

  dd <- vector(mode = "numeric", length = loop_size)
  for (jj in 1:loop_size) {
    reg_object_boot <- NULL
    while (is.null(reg_object_boot)) {
      boot_samp <- sample(1:n_nons, size = n_nons, replace = TRUE)
      # boot_samp <- sample(1:n_rand, size = n_rand, replace = TRUE)
      y_nons_b <- y[boot_samp]

      reg_object_boot <- switch(pmm_reg_engine,
        "glm" = stats::glm(
          formula = model_obj$model$glm_object$formula,
          data = model_obj$model$glm_object$data[boot_samp, , drop = FALSE],
          # weights = weights,
          family = model_obj$model$glm_object$family,
          start = model_obj$model$glm_object$coefficients
        ),
        "loess" = stats::loess(
          formula = model_obj$model$glm_object$formula,
          data = model_obj$model$glm_object$data[boot_samp, , drop = FALSE],
          span = .2,
          control = stats::loess.control(surface = "direct")
        )
      )
      XX <- predict(
        reg_object_boot,
        newdata = svydesign$variables,
        type = "response"
      )
      # XX <- reg_object_boot$family$mu.eta(X_rand %*% reg_object_boot$coefficients)

      if (any(!is.finite(XX))) {
        reg_object_boot <- NULL
      }
    }

    YY <- switch(predictive_match,
      {
        nonprobMI_nn(
          data = y_nons_b,
          query = XX,
          k = k,
          searchtype = "standard",
          treetype = "kd"
        )
      },
      {
        nonprobMI_nn(
          data = predict(
            reg_object_boot,
            newdata = model_obj$model$glm_object$data[boot_samp, , drop = FALSE],
            type = "response"
          ),
          query = XX,
          k = k,
          searchtype = "standard",
          treetype = "kd"
        )
      }
    )

    dd[jj] <- weighted.mean(
      apply(YY$nn.idx, 1, FUN = \(x) mean(y_nons_b[x])),
      weights_rand
    )
  }
  var(dd)
}
