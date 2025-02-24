#' @useDynLib nonprobsvy
#' @import Rcpp
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats binomial
#' @importFrom stats terms
#' @importFrom stats reformulate
#' @importFrom survey svytotal
#' @importFrom ncvreg cv.ncvreg
#' @importFrom MASS ginv
#' @importFrom Rcpp evalCpp
#' @importFrom formula.tools rhs
#' @importFrom formula.tools lhs
#' @noRd
nonprob_dr <- function(selection,
                       outcome,
                       data,
                       svydesign,
                       pop_totals,
                       pop_means,
                       pop_size,
                       method_selection,
                       method_outcome,
                       family_outcome,
                       subset,
                       strata,
                       weights,
                       na_action,
                       control_selection,
                       control_outcome,
                       control_inference,
                       start_outcome,
                       start_selection,
                       verbose,
                       se,
                       pop_size_fixed,
                       ...) {

  # important parameters:
  #   bias_correction = FALSE,
  #   bias_inf = c("union", "div"),

  ## setting for the IPW
  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  estimation_method <- est_method_ipw(control_selection$est_method)

  ## multiple outcomes
  outcomes <- make_outcomes(outcome)
  output <- list()
  outcome_list <- list()
  ys <- list()

  confidence_interval <- list()
  SE_values <- list()


  ###
  if (control_inference$vars_selection & control_inference$vars_combine) {

    ## estimate the mi
    if (verbose) {
      cat("MI variable selection in progress...\n")
    }
    results_mi <- nonprob_mi(outcome = outcome,
                             data = data,
                             svydesign = svydesign,
                             pop_totals = pop_totals,
                             pop_means = pop_means,
                             pop_size = pop_size,
                             method_outcome = method_outcome,
                             family_outcome = family_outcome,
                             subset = subset,
                             strata = strata,
                             weights = weights,
                             na_action = na_action,
                             control_outcome = control_outcome,
                             control_inference = control_inference,
                             start_outcome = start_outcome,
                             verbose = verbose,
                             se = FALSE,
                             pop_size_fixed=pop_size_fixed)

    if (verbose) {
      cat("IPW variable selection in progress...\n")
    }
    results_ipw <- nonprob_ipw(selection = selection,
                               target = reformulate(outcomes[[1]]),
                               data = data,
                               svydesign = svydesign,
                               pop_totals = pop_totals,
                               pop_means = pop_means,
                               pop_size = pop_size,
                               method_selection = method_selection,
                               subset = subset,
                               strata = strata,
                               weights = weights,
                               na_action = na_action,
                               control_selection = control_selection,
                               control_inference = control_inference,
                               start_selection = start_selection,
                               verbose = verbose,
                               se = FALSE,
                               pop_size_fixed = pop_size_fixed)

    ## doubly robust estimator
    mu_hat <- mu_hatDR(y_hat = results_mi$output$mean,
                       y_resid = do.call("cbind", results_mi$ys_resid),
                       weights = weights,
                       weights_nons = results_ipw$ipw_weights,
                       N_nons = sum(results_ipw$ipw_weights))

    ipw_coefs_sel <- names(results_ipw$selection$coefficients)
    mi_coefs_sel <- lapply(results_mi$outcome, coef)
    dr_coefs_sel <- lapply(mi_coefs_sel, function(x) {
      mi_cols <- names(x[abs(x)>0])
      combined <- sort(base::union(ipw_coefs_sel, mi_cols))
      combined[!grepl("Intercept", combined)]
    })

    if (verbose) {
      cat("IPW vars selected:", ipw_coefs_sel, "\n")
      cat("MI vars selected:\n")
      print(lapply(mi_coefs_sel, function(x) names(x[abs(x)>0])))
      cat("DR combined vars:\n")
      print(dr_coefs_sel)
    }


    selection_vars <- all.vars(formula.tools::rhs(outcome))
    outcome_vars <- all.vars(formula.tools::rhs(selection))
    target_vars <- all.vars(formula.tools::lhs(outcome))
    combined_vars <- reformulate(union(selection_vars, outcome_vars))

    y_nons <- subset(data, select=target_vars)
    X_nons <- model.matrix(combined_vars, data)
    X_rand <- model.matrix(combined_vars, svydesign$variables) ## if design is present

    X_nons <- cbind(y_nons, X_nons[, !grepl("Intercept", colnames(X_nons)), drop=FALSE])
    X_rand <- X_rand[, !grepl("Intercept", colnames(X_rand)), drop=FALSE]
    svydesign_ <- svydesign
    svydesign_$variables <- cbind(svydesign_$variables, X_rand)

    results_mi_combined <- results_ipw_combined <- list()
    mu_hat <- numeric()

    control_inference_ <- control_inference
    control_inference_$vars_selection <- FALSE

    for (o in outcomes$f) {
      ## this is not saved in the output list
      results_ipw_combined[[o]] <- nonprob_ipw(data = X_nons,
                                               target = reformulate(o),
                                               selection =  reformulate(dr_coefs_sel[[o]]),
                                               svydesign = svydesign_,
                                               pop_totals = pop_totals,
                                               pop_means = pop_means,
                                               pop_size = pop_size,
                                               method_selection = method_selection,
                                               subset = subset,
                                               strata = strata,
                                               weights = weights,
                                               na_action = na_action,
                                               control_selection = control_selection,
                                               control_inference = control_inference_,
                                               start_selection = start_selection,
                                               verbose = verbose,
                                               se = FALSE,
                                               pop_size_fixed = pop_size_fixed)
      ## estimate the mi
      results_mi_combined[[o]] <- nonprob_mi(outcome = as.formula(paste0(o, reformulate(dr_coefs_sel[[o]]))),
                                             data = X_nons,
                                             svydesign = svydesign_,
                                             pop_totals = pop_totals,
                                             pop_means = pop_means,
                                             pop_size = pop_size,
                                             method_outcome = method_outcome,
                                             family_outcome = family_outcome,
                                             subset = subset,
                                             strata = strata,
                                             weights = weights,
                                             na_action = na_action,
                                             control_outcome = control_outcome,
                                             control_inference = control_inference_,
                                             start_outcome = start_outcome,
                                             verbose = verbose,
                                             se = FALSE,
                                             pop_size_fixed=pop_size_fixed)

      ## estimate in loop
      mu_hat[[o]] <- mu_hatDR(y_hat = results_mi_combined[[o]]$output$mean,
                              y_resid = do.call("cbind", results_mi_combined[[o]]$ys_resid),
                              weights = weights,
                              weights_nons = results_ipw_combined[[o]]$ipw_weights,
                              N_nons = sum(results_ipw_combined[[o]]$ipw_weights))
    }

    if (se) {
      for (o in outcomes$f) {
        if (control_inference$var_method == "analytic") {
            b_var <- method$b_vec_dr(X = results_ipw_combined[[o]]$X[results_ipw_combined[[o]]$R == 1, ],
                                     ps = results_ipw_combined[[o]]$ps_scores[results_ipw_combined[[o]]$R == 1],
                                     psd = as.numeric(results_ipw_combined[[o]]$selection$selection_model$ps_nons_der),
                                     y = results_mi_combined[[o]]$y[[o]],
                                     hess = results_ipw_combined[[o]]$selection$selection_model$hess,
                                     eta = as.numeric(results_ipw_combined[[o]]$X[results_ipw_combined[[o]]$R == 1, ] %*% as.matrix(results_ipw_combined[[o]]$selection$coefficients)),
                                     h_n = 1 / pop_size * sum(results_mi_combined[[o]]$y[[o]] - results_mi_combined[[o]]$ys_nons_pred[[1]]),
                                     y_pred = results_mi_combined[[o]]$ys_nons_pred[[1]],
                                     weights = weights,
                                     verbose = verbose)

            var_nonprob <- estimation_method$make_var_nonprob(
              ps = results_ipw_combined[[o]]$ps_scores[results_ipw_combined[[o]]$R == 1],
              psd = as.numeric(results_ipw_combined[[o]]$selection$selection_model$ps_nons_der),
              y = results_mi_combined[[o]]$y[[o]],
              y_pred = results_mi_combined[[o]]$ys_nons_pred[[o]],
              h_n = 1 / pop_size * sum(results_mi_combined[[o]]$y[[o]] - results_mi_combined[[o]]$ys_nons_pred[[1]]),
              X = results_ipw_combined[[o]]$X[results_ipw_combined[[o]]$R == 1, ],
              b = b_var,
              N = pop_size,
              gee_h_fun = control_selection$gee_h_fun,
              method_selection = method_selection,
              weights = weights,
              pop_totals = pop_totals
            )

            if (is.null(pop_totals)) {
              t_comp <- estimation_method$make_t_comp(
                X = results_ipw_combined[[o]]$X[results_ipw_combined[[o]]$R == 0, ],
                ps = as.numeric(results_ipw_combined[[o]]$selection$selection_model$est_ps_rand),
                psd = as.numeric(results_ipw_combined[[o]]$selection$selection_model$est_ps_rand_der),
                b = b_var,
                gee_h_fun = control_selection$gee_h_fun,
                y_rand = results_mi_combined[[o]]$ys_rand_pred[[1]],
                y_nons = results_mi_combined[[o]]$ys_nons_pred[[1]],
                N = pop_size,
                method_selection = method_selection,
                weights = weights
              )

              svydesign__ <- stats::update(svydesign_, t_comp = t_comp)
              svydesign_mean <- survey::svymean(~t_comp, svydesign__)
              var_prob <- as.vector(attr(svydesign_mean, "var"))

            } else {
              var_prob <- 0
            }

        var_total <- var_nonprob + var_prob
        se_nonprob <- sqrt(var_nonprob)
        se_prob <- sqrt(var_prob)
        SE_values[[o]] <- data.frame(prob = se_prob, nonprob = se_nonprob)
        z <- stats::qnorm(1 - control_inference$alpha / 2)
        # confidence interval based on the normal approximation
        confidence_interval[[o]] <- data.frame(lower_bound = mu_hat - z * sqrt(var_total),
                                               upper_bound = mu_hat + z * sqrt(var_total))
        output[[o]] <- data.frame(mean = mu_hat[o], SE = sqrt(var_total))
        }
      }
        ## this will not work but needs to be updated
        if (control_inference$var_method == "bootstrap") {

        boot_obj <- boot_dr(selection = selection,
                            outcome = outcome,
                            target = reformulate(outcomes[[1]]),
                            data = data,
                            svydesign = svydesign,
                            pop_totals = pop_totals,
                            pop_means = pop_means,
                            pop_size = pop_size,
                            method_selection = method_selection,
                            method_outcome = method_outcome,
                            family_outcome = family_outcome,
                            subset = subset,
                            strata = strata,
                            weights = weights,
                            na_action = na_action,
                            control_selection = control_selection,
                            control_outcome = control_outcome,
                            control_inference = control_inference,
                            start_outcome = start_outcome,
                            start_selection = start_selection,
                            verbose = verbose,
                            pop_size_fixed = pop_size_fixed)

        var_total <- apply(boot_obj, 2, var)
        SE_values <- replicate(NROW(outcomes[[1]]), data.frame(nonprob = NA, prob = NA), simplify = F)
        SE <- sqrt(var_total)
        output <- list(data.frame(mean = mu_hat, SE = SE))
        alpha <- control_inference$alpha
        z <- stats::qnorm(1 - alpha / 2)
        # confidence interval based on the normal approximation
        confidence_interval <- list(data.frame(lower_bound = mu_hat - z * SE,
                                               upper_bound = mu_hat + z * SE))
        }
      } else {
      for (o in 1:outcomes$l) {
        confidence_interval[[o]] <- data.frame(lower_bound = NA, upper_bound = NA)
        SE_values[[o]] <- data.frame(nonprob = NA, prob = NA)
        output[[o]] <- data.frame(mean = mu_hat[o], SE = NA)
      }
    }
  } else {
    ## without variable selection
    results_mi <- nonprob_mi(outcome = outcome,
                             data = data,
                             svydesign = svydesign,
                             pop_totals = pop_totals,
                             pop_means = pop_means,
                             pop_size = pop_size,
                             method_outcome = method_outcome,
                             family_outcome = family_outcome,
                             subset = subset,
                             strata = strata,
                             weights = weights,
                             na_action = na_action,
                             control_outcome = control_outcome,
                             control_inference = control_inference,
                             start_outcome = start_outcome,
                             verbose = verbose,
                             se = FALSE,
                             pop_size_fixed=pop_size_fixed)

    results_ipw <- nonprob_ipw(selection = selection,
                               target = reformulate(outcomes[[1]]),
                               data = data,
                               svydesign = svydesign,
                               pop_totals = pop_totals,
                               pop_means = pop_means,
                               pop_size = pop_size,
                               method_selection = method_selection,
                               subset = subset,
                               strata = strata,
                               weights = weights,
                               na_action = na_action,
                               control_selection = control_selection,
                               control_inference = control_inference,
                               start_selection = start_selection,
                               verbose = verbose,
                               se = FALSE,
                               pop_size_fixed = pop_size_fixed)

    ## doubly robust estimator
    mu_hat <- mu_hatDR(y_hat = results_mi$output$mean,
                       y_resid = do.call("cbind", results_mi$ys_resid),
                       weights = weights,
                       weights_nons = results_ipw$ipw_weights,
                       N_nons = sum(results_ipw$ipw_weights))

    if (se) {
      if (control_inference$var_method == "analytic") {
        for (o in 1:outcomes$l) {
          b_var <- method$b_vec_dr(X = results_ipw$X[results_ipw$R == 1, ],
                                   ps = results_ipw$ps_scores[results_ipw$R == 1],
                                   psd = as.numeric(results_ipw$selection$selection_model$ps_nons_der),
                                   y = results_mi$y[[o]],
                                   hess = results_ipw$selection$selection_model$hess,
                                   eta = as.numeric(results_ipw$X[results_ipw$R == 1, ] %*% as.matrix(results_ipw$selection$coefficients)),
                                   h_n = 1 / pop_size * sum(results_mi$y[[o]] - results_mi$ys_nons_pred[[o]]),
                                   y_pred = results_mi$ys_nons_pred[[o]],
                                   weights = weights,
                                   verbose = verbose)

          var_nonprob <- estimation_method$make_var_nonprob(
            ps = results_ipw$ps_scores[results_ipw$R == 1],
            psd = as.numeric(results_ipw$selection$selection_model$ps_nons_der),
            y = results_mi$y[[o]],
            y_pred = results_mi$ys_nons_pred[[o]],
            h_n = 1 / pop_size * sum(results_mi$y[[o]] - results_mi$ys_nons_pred[[o]]),
            X = results_ipw$X[results_ipw$R == 1, ],
            b = b_var,
            N = pop_size,
            gee_h_fun = control_selection$gee_h_fun,
            method_selection = method_selection,
            weights = weights,
            pop_totals = pop_totals
          )

          if (is.null(pop_totals)) {
            t_comp <- estimation_method$make_t_comp(
              X = results_ipw$X[results_ipw$R == 0, ],
              ps = as.numeric(results_ipw$selection$selection_model$est_ps_rand),
              psd = as.numeric(results_ipw$selection$selection_model$est_ps_rand_der),
              b = b_var,
              gee_h_fun = control_selection$gee_h_fun,
              y_rand = results_mi$ys_rand_pred[[o]],
              y_nons = results_mi$ys_nons_pred[[o]],
              N = pop_size,
              method_selection = method_selection,
              weights = weights
            )

            svydesign_ <- stats::update(svydesign, t_comp = t_comp)
            svydesign_mean <- survey::svymean(~t_comp, svydesign)
            var_prob <- as.vector(attr(svydesign_mean, "var"))

          } else {
            var_prob <- 0
          }

          var_total <- var_nonprob + var_prob
          se_nonprob <- sqrt(var_nonprob)
          se_prob <- sqrt(var_prob)
          SE_values[[o]] <- data.frame(prob = se_prob, nonprob = se_nonprob)
          z <- stats::qnorm(1 - control_inference$alpha / 2)
          # confidence interval based on the normal approximation
          confidence_interval[[o]] <- data.frame(lower_bound = mu_hat - z * sqrt(var_total),
                                                 upper_bound = mu_hat + z * sqrt(var_total))
          output[[o]] <- data.frame(mean = mu_hat[o], SE = sqrt(var_total))
        }
      }

      if (control_inference$var_method == "bootstrap") {

        boot_obj <- boot_dr(selection = selection,
                            outcome = outcome,
                            target = reformulate(outcomes[[1]]),
                            data = data,
                            svydesign = svydesign,
                            pop_totals = pop_totals,
                            pop_means = pop_means,
                            pop_size = pop_size,
                            method_selection = method_selection,
                            method_outcome = method_outcome,
                            family_outcome = family_outcome,
                            subset = subset,
                            strata = strata,
                            weights = weights,
                            na_action = na_action,
                            control_selection = control_selection,
                            control_outcome = control_outcome,
                            control_inference = control_inference,
                            start_outcome = start_outcome,
                            start_selection = start_selection,
                            verbose = verbose,
                            pop_size_fixed = pop_size_fixed)

        var_total <- apply(boot_obj, 2, var)
        SE_values <- replicate(NROW(outcomes[[1]]), data.frame(nonprob = NA, prob = NA), simplify = F)
        SE <- sqrt(var_total)
        output <- list(data.frame(mean = mu_hat, SE = SE))
        alpha <- control_inference$alpha
        z <- stats::qnorm(1 - alpha / 2)
        # confidence interval based on the normal approximation
        confidence_interval <- list(data.frame(lower_bound = mu_hat - z * SE,
                                               upper_bound = mu_hat + z * SE))
      }

    } else {
      for (o in 1:outcomes$l) {
        confidence_interval[[o]] <- data.frame(lower_bound = NA, upper_bound = NA)
        SE_values[[o]] <- data.frame(nonprob = NA, prob = NA)
        output[[o]] <- data.frame(mean = mu_hat[o], SE = NA)
      }
    }


  }
  boot_sample <- if (control_inference$var_method == "bootstrap" & control_inference$keep_boot) {
    boot_obj
  } else {
    NULL
  }

  if (!is.null(boot_sample) & is.matrix(boot_sample)) colnames(boot_sample) <- names(ys)

  output <- do.call(rbind, output)
  confidence_interval <- do.call(rbind, confidence_interval)
  SE_values <- do.call(rbind, SE_values)
  rownames(output) <- rownames(confidence_interval) <- rownames(SE_values) <- outcomes$f

  structure(
    list(
      data = data,
      X = results_mi$X,
      y = results_mi$y,
      R = results_ipw$R,
      ps_scores = results_ipw$prop_scores,
      case_weights = results_ipw$case_weights,
      ipw_weights = results_ipw$ipw_weights,
      control = list(
        control_selection = control_selection,
        control_outcome = control_outcome,
        control_inference = control_inference
      ),
      output = output,
      SE = SE_values,
      confidence_interval = confidence_interval,
      nonprob_size = results_mi$nonprob_size,
      prob_size = results_mi$prob_size,
      pop_size = pop_size,
      pop_size_fixed = pop_size_fixed,
      pop_totals = pop_totals,
      pop_means = pop_means,
      outcome = results_mi$outcome,
      selection = results_ipw$selection,
      boot_sample = boot_sample,
      svydesign = if (is.null(pop_totals)) svydesign else NULL,
      ys_rand_pred = results_mi$ys_rand_pred,
      ys_nons_pred = results_mi$ys_nons_pred,
      ys_resid = results_mi$ys_resid
    ),
    class = "nonprob"
  )
}
