#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats weighted.mean
#' @importFrom RANN nn2
#' @importFrom stats terms
#' @rdname main_doc

nonprobMI <- function(outcome,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_outcome,
                      control_inference = controlInf(var_method = "analytic"),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  outcomes <- ff(outcome)
  output <- list()
  confidence_interval <- list()
  SE_values <- list()
  num_boot <- control_inference$num_boot
  for (k in 1:outcomes$l) {
    if (is.null(pop_totals) && !is.null(svydesign)) {
      outcome <- outcomes$outcome[[k]]

      # model for outcome formula
      OutcomeModel <- model_frame(formula = outcome, data = data, svydesign = svydesign)
      X_nons <- OutcomeModel$X_nons
      X_rand <- OutcomeModel$X_rand
      nons_names <- OutcomeModel$nons_names
      y_nons <- OutcomeModel$y_nons

      R_nons <- rep(1, nrow(X_nons))
      R_rand <- rep(0, nrow(X_rand))
      R <- c(R_nons, R_rand)

      loc_nons <- which(R == 1)
      loc_rand <- which(R == 0)

      n_nons <- nrow(X_nons)
      n_rand <- nrow(X_rand)
      X <- rbind(X_nons, X_rand)

      ps_rand <- svydesign$prob
      weights_rand <- 1/ps_rand
      N_est_rand <- sum(weights_rand)

      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      ## estimation
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(outcome = outcome,
                                 data = data,
                                 weights = weights,
                                 family_outcome = family_outcome,
                                 X_nons = X_nons,
                                 y_nons = y_nons,
                                 X_rand = X_rand,
                                 control = control_outcome,
                                 n_nons = n_nons,
                                 n_rand = n_rand,
                                 model_frame = OutcomeModel$model_frame_rand,
                                 vars_selection = control_inference$vars_selection,
                                 pop_totals = pop_totals)
      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      model_out <- model_obj$model
      parameters <- model_obj$parameters

      # updating probability sample by adding y_hat variable
      svydesign <- stats::update(svydesign,
                                 y_hat_MI = y_rand_pred)
      mu_hat <- weighted.mean(y_rand_pred, w = weights_rand)
    } else if ((!is.null(pop_totals) || !is.null(pop_means)) && is.null(svydesign)) {

      if (!is.null(pop_totals)) pop_size <- N_est_rand <- pop_totals[1]
      if (!is.null(pop_means)) { # TO consider
        if (!is.null(pop_size)) {
          pop_totals <- c(pop_size, pop_size * pop_means)
          names(pop_totals) <- c("(Intercept)", names(pop_means))
        } else {
          stop("pop_size must be defined when estimating with pop_means.")
        }
      }
      Model <- model_frame(formula = outcome, data = data, pop_totals = pop_totals)

      X_nons <- Model$X_nons
      X_rand <- Model$X_rand
      y_nons <- Model$y_nons
      n_rand <- 0
      weights_rand <- NULL
      n_nons <- nrow(X_nons)

      method_outcome_nonprobsvy <- paste(method_outcome, "_nonprobsvy", sep = "")
      MethodOutcome <- get(method_outcome_nonprobsvy, mode = "function", envir = parent.frame())
      model_obj <- MethodOutcome(outcome = outcome,
                                 data = data,
                                 weights = weights,
                                 family_outcome = family_outcome,
                                 X_nons = X_nons,
                                 y_nons = y_nons,
                                 X_rand = X_rand,
                                 control = control_outcome,
                                 n_nons = n_nons,
                                 n_rand = n_rand,
                                 model_frame = OutcomeModel$model_frame_rand,
                                 vars_selection = control_inference$vars_selection,
                                 pop_totals = pop_totals)

      y_rand_pred <- model_obj$y_rand_pred
      y_nons_pred <- model_obj$y_nons_pred
      model_out <- model_obj$model
      parameters <- model_obj$parameters

      mu_hat <- ifelse(method_outcome == "glm", as.vector(y_rand_pred/N_est_rand), y_rand_pred)
    } else {
      stop("Please, provide svydesign object or pop_totals/pop_means.")
    }

    # design based variance estimation based on approximations of the second-order inclusion probabilities
    if (control_inference$var_method == "analytic") { # consider move variance implementation to internals
      var_obj <- internal_varMI(svydesign = svydesign,
                                X_nons = X_nons,
                                X_rand = X_rand,
                                y = y_nons,
                                y_pred = y_nons_pred,
                                weights_rand = weights_rand,
                                method = method_outcome,
                                n_rand = n_rand,
                                n_nons = n_nons,
                                N = N_est_rand,
                                family = family_outcome,
                                parameters = parameters,
                                pop_totals = pop_totals)

      var_nonprob <- var_obj$var_nonprob
      var_prob <- var_obj$var_prob

      se_nonprob <- sqrt(var_nonprob)
      se_prob <- sqrt(var_prob)
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(prob = se_prob, nonprob = se_nonprob))))
      # variance
      var <- var_nonprob + var_prob

    } else if (control_inference$var_method == "bootstrap") { # TODO for pop_totals
      # bootstrap variance
      if (control_inference$cores > 1) {
        var <- bootMI_multicore(X_rand,
                                X_nons,
                                weights,
                                y_nons,
                                family_outcome,
                                num_boot = num_boot,
                                weights_rand,
                                mu_hat,
                                svydesign,
                                rep_type = control_inference$rep_type,
                                method = method_outcome,
                                control = control_outcome,
                                pop_totals = pop_totals,
                                cores = control_inference$cores,
                                verbose = verbose)
      } else {
        var <- bootMI(X_rand,
                      X_nons,
                      weights,
                      y_nons,
                      family_outcome,
                      num_boot = num_boot,
                      weights_rand,
                      mu_hat,
                      svydesign,
                      rep_type = control_inference$rep_type,
                      method = method_outcome,
                      control = control_outcome,
                      pop_totals = pop_totals,
                      verbose = verbose)
      }
      SE_values[[k]] <- data.frame(t(data.frame("SE" = c(nonprob = "no division into nonprobability", prob = "probability sample in case of bootstrap variance"))))
    }

    X <- rbind(X_nons, X_rand) # joint model matrix
    #if (is.null(pop_size))
    pop_size <- N_est_rand # estimated pop_size
    se <- sqrt(var)

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    confidence_interval[[k]] <- data.frame(t(data.frame("normal" = c(lower_bound = mu_hat - z * se,
                                                                upper_bound = mu_hat + z * se
    ))))

    output[[k]] <- data.frame(t(data.frame(result = c(mean = mu_hat, SE = se))))
  }
  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f

  structure(
    list(X = X,
         control = list(control_outcome = control_outcome,
                        control_inference = control_inference,
                        method_outcome = method_outcome),
         output = output,
         SE_values = SE_values,
         confidence_interval = confidence_interval,
         parameters = parameters,
         nonprob_size = n_nons,
         prob_size = n_rand,
         pop_size = pop_size,
         outcome = model_out
         ),
    class = c("nonprobsvy", "nonprobsvy_mi"))

  }
