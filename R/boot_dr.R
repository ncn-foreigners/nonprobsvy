#' Bootstrap for the DR estimator
#' @noRd
boot_dr <- function(selection,
                    outcome,
                    target,
                    data,
                    svydesign,
                    pop_totals,
                    pop_means,
                    pop_size,
                    method_selection,
                    method_outcome,
                    family_outcome,
                    strata,
                    case_weights,
                    na_action,
                    control_selection,
                    control_outcome,
                    control_inference,
                    start_outcome,
                    start_selection,
                    verbose,
                    pop_size_fixed) {

  # Initialize objects to store results
  num_boot <- control_inference$num_boot
  bias_corr <- control_inference$bias_correction
  vars_combine <- control_inference$vars_combine

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  ## add bootstrap with variables combination

  if (!is.null(svydesign)) {
    svydesign_rep <- survey::as.svrepdesign(svydesign,
                                            type = control_inference$rep_type,
                                            replicates = num_boot)
    rep_weights <- svydesign_rep$repweights$weights
  }

  # Single core processing
  if (control_inference$cores == 1) {

    boot_obj <- matrix(0, ncol = length(all.vars(target)), nrow = num_boot)

    if (verbose) {
      message("Single core bootstrap in progress...")
      pb_boot <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
    }

    b <- 1

    if (!is.null(svydesign)) {

      # Bootstrap for probability and non-probability samples
      while (b <= num_boot) {

        ## probability part
        strap_rand_svy <- which(rep_weights[, b] != 0)
        weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
        pop_size_strap <- sum(weights_rand_strap_svy)

        data_prob <- svydesign$variables[strap_rand_svy, ]
        data_prob$weight <- weights_rand_strap_svy[strap_rand_svy]

        svyd_call <- svydesign$call
        svyd_call$data <- as.name("data_prob")
        svydesign_b <- eval(svyd_call)

        strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / case_weights)

        tryCatch(
          {
            results_ipw_b <- nonprob_ipw(selection = selection,
                                         target = target,
                                         data = data[strap_nons, ],
                                         svydesign = svydesign_b,
                                         pop_totals = NULL,
                                         pop_means = NULL,
                                         pop_size = NULL,
                                         method_selection = method_selection,
                                         strata = strata,
                                         case_weights = case_weights[strap_nons],
                                         na_action = na_action,
                                         control_selection = control_selection,
                                         control_inference = control_inference,
                                         start_selection = start_selection,
                                         verbose = FALSE,
                                         se = FALSE,
                                         pop_size_fixed=pop_size_fixed)
            ## estimate the mi
            results_mi_b <- nonprob_mi(outcome = outcome,
                                       data = data[strap_nons, ],
                                       svydesign = svydesign_b,
                                       pop_totals = NULL,
                                       pop_means = NULL,
                                       pop_size = NULL,
                                       method_outcome = method_outcome,
                                       family_outcome = family_outcome,
                                       strata = strata,
                                       case_weights = case_weights[strap_nons],
                                       na_action = na_action,
                                       control_outcome = control_outcome,
                                       control_inference = control_inference,
                                       start_outcome = start_outcome,
                                       verbose = FALSE,
                                       se = FALSE,
                                       pop_size_fixed=pop_size_fixed)

            ## combination of variables after bootstrap
            if (!vars_combine) {
              boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                        y_resid = do.call("cbind", results_mi_b$ys_resid),
                                        weights = case_weights[strap_nons],
                                        weights_nons = results_ipw_b$ipw_weights,
                                        N_nons = sum(results_ipw_b$ipw_weights))
            } else {
              if (verbose) message("\nCombining variables...")

              ipw_coefs_sel <- names(results_ipw_b$selection$coefficients)
              mi_coefs_sel <- lapply(results_mi_b$outcome, coef)
              dr_coefs_sel <- lapply(mi_coefs_sel, function(x) {
                mi_cols <- names(x[abs(x)>0])
                combined <- sort(base::union(ipw_coefs_sel, mi_cols))
                combined[!grepl("Intercept", combined)]
              })

              ## combining variables for selection
              selection_vars <- all.vars(formula.tools::rhs(outcome))
              outcome_vars <- all.vars(formula.tools::rhs(selection))
              target_vars <- all.vars(formula.tools::lhs(outcome))
              combined_vars <- reformulate(union(selection_vars, outcome_vars))

              y_nons <- subset(data[strap_nons, ], select=target_vars)
              X_nons <- model.matrix(combined_vars, data[strap_nons, ])
              X_rand <- model.matrix(combined_vars, svydesign_b$variables) ## if design is present
              X_all <- rbind(X_rand, X_nons)

              X_nons <- cbind(y_nons, X_nons[, !grepl("Intercept", colnames(X_nons)), drop=FALSE])
              X_rand <- X_rand[, !grepl("Intercept", colnames(X_rand)), drop=FALSE]
              svydesign_b_ <- svydesign_b
              svydesign_b_$variables <- cbind(svydesign_b_$variables, X_rand)

              if (bias_corr) {

                if (verbose) message("\nBias correction...")

                ## consider different start
                par0 <- numeric(NCOL(X_all)*2)
                names(par0) <- rep(colnames(X_all), times = 2)

                bias_corr_result_b <- nleqslv::nleqslv(
                  x = par0,
                  fn = u_theta_beta_dr,
                  method = control_selection$nleqslv_method,
                  global = control_selection$nleqslv_global,
                  xscalm = control_selection$nleqslv_xscalm,
                  jacobian = TRUE,
                  control = list(
                    scalex = rep(1, length(par0)),
                    maxit = control_selection$maxit
                  ),
                  R = results_ipw_b$R,
                  X = X_all,
                  y = c(rep(0, sum(results_ipw_b$R==0)), y_nons[, 1]),
                  weights = c(weights(svydesign_b_), results_ipw_b$case_weights),#c(weights(svydesign_), weights),
                  method_selection = method_selection,
                  family_outcome = family_outcome
                )

                theta_hat <- bias_corr_result_b$x[1:NCOL(X_all)]
                beta_hat <- bias_corr_result_b$x[(NCOL(X_all) + 1):(2 * NCOL(X_all))]

                bias_corr_ps <- method$make_link_inv(unname(drop(X_all %*% theta_hat)))
                bias_corr_ipw_weights <- 1/bias_corr_ps[results_ipw_b$R == 1]
                bias_corr_mu_rand_pred <- as.vector(get(family_outcome)()$linkinv(X_all[results_ipw_b$R == 0, ] %*% beta_hat))
                bias_corr_mu_nons_pred <- as.vector(get(family_outcome)()$linkinv(X_all[results_ipw_b$R == 1, ] %*% beta_hat))
                bias_corr_mu_resid <- bias_corr_mu_nons_pred - y_nons

                boot_obj[b, ] <- mu_hatDR(y_hat = weighted.mean(bias_corr_mu_rand_pred, weights(svydesign_b_)),
                                          y_resid = as.matrix(bias_corr_mu_resid),
                                          weights = results_ipw_b$case_weights,
                                          weights_nons = bias_corr_ipw_weights,
                                          N_nons = sum(bias_corr_ipw_weights))
              } else {
                results_ipw_b_combined <- nonprob_ipw(data = as.data.frame(X_nons),
                                                      target = outcome,
                                                      selection =  reformulate(dr_coefs_sel[[1]]),
                                                      svydesign = svydesign_b_,
                                                      pop_totals = NULL,
                                                      pop_means = NULL,
                                                      pop_size = NULL,
                                                      method_selection = method_selection,
                                                      strata = strata,
                                                      case_weights = case_weights[strap_nons],
                                                      na_action = na_action,
                                                      control_selection = control_selection,
                                                      control_inference = control_inference,
                                                      start_selection = start_selection,
                                                      verbose = FALSE,
                                                      se = FALSE,
                                                      pop_size_fixed = pop_size_fixed)
                ## estimate the mi
                results_mi_b_combined <- nonprob_mi(outcome = as.formula(paste0(target_vars, reformulate(dr_coefs_sel[[1]]))),
                                                    data = as.data.frame(X_nons),
                                                    svydesign = svydesign_b_,
                                                    pop_totals = NULL,
                                                    pop_means = NULL,
                                                    pop_size = NULL,
                                                    method_outcome = method_outcome,
                                                    family_outcome = family_outcome,
                                                    strata = strata,
                                                    case_weights = case_weights[strap_nons],
                                                    na_action = na_action,
                                                    control_outcome = control_outcome,
                                                    control_inference = control_inference,
                                                    start_outcome = start_outcome,
                                                    verbose = FALSE,
                                                    se = FALSE,
                                                    pop_size_fixed=pop_size_fixed)

                boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b_combined$output$mean,
                                          y_resid = do.call("cbind", results_mi_b_combined$ys_resid),
                                          weights = case_weights[strap_nons],
                                          weights_nons = results_ipw_b_combined$ipw_weights,
                                          N_nons = sum(results_ipw_b_combined$ipw_weights))

              }
            }
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", b, " iteration: ", e$message, sep = "")
              message(info)
            }
          }
        )
        if (verbose) {
          utils::setTxtProgressBar(pb_boot, b)
        }
        b <- b + 1
      }
    } else {
      # Bootstrap for non-probability samples only
      while (b <= num_boot) {

        strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / case_weights)

        tryCatch(
          {
            # Non-probability part
            results_ipw_b <- nonprob_ipw(selection = selection,
                                         target = target,
                                         data = data[strap_nons, ],
                                         svydesign = svydesign,
                                         pop_totals = pop_totals,
                                         pop_means = pop_means,
                                         pop_size = pop_size,
                                         method_selection = method_selection,
                                         strata = strata,
                                         case_weights = case_weights[strap_nons],
                                         na_action = na_action,
                                         control_selection = control_selection,
                                         control_inference = control_inference,
                                         start_selection = start_selection,
                                         verbose = FALSE,
                                         se = FALSE,
                                         pop_size_fixed=pop_size_fixed)
            ## estimate the mi
            results_mi_b <- nonprob_mi(outcome = outcome,
                                       data = data[strap_nons, ],
                                       svydesign = svydesign,
                                       pop_totals = pop_totals,
                                       pop_means = pop_means,
                                       pop_size = pop_size,
                                       method_outcome = method_outcome,
                                       family_outcome = family_outcome,
                                       strata = strata,
                                       case_weights = case_weights[strap_nons],
                                       na_action = na_action,
                                       control_outcome = control_outcome,
                                       control_inference = control_inference,
                                       start_outcome = start_outcome,
                                       verbose = FALSE,
                                       se = FALSE,
                                       pop_size_fixed=pop_size_fixed)

            ## combination of variables after bootstrap
            if (!vars_combine) {
              boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                        y_resid = do.call("cbind", results_mi_b$ys_resid),
                                        weights = case_weights[strap_nons],
                                        weights_nons = results_ipw_b$ipw_weights,
                                        N_nons = sum(results_ipw_b$ipw_weights))
            } else {

              if (verbose) message("\nCombining variables...")

              ipw_coefs_sel <- names(results_ipw_b$selection$coefficients)
              mi_coefs_sel <- lapply(results_mi_b$outcome, coef)
              dr_coefs_sel <- lapply(mi_coefs_sel, function(x) {
                mi_cols <- names(x[abs(x)>0])
                combined <- sort(base::union(ipw_coefs_sel, mi_cols))
                combined[!grepl("Intercept", combined)]
              })

              ## combining variables for selection
              selection_vars <- all.vars(formula.tools::rhs(outcome))
              outcome_vars <- all.vars(formula.tools::rhs(selection))
              target_vars <- all.vars(formula.tools::lhs(outcome))
              combined_vars <- reformulate(union(selection_vars, outcome_vars))

              y_nons <- subset(data[strap_nons, ], select=target_vars)
              X_nons <- model.matrix(combined_vars, data[strap_nons, ])
              pop_totals_ <- pop_totals[colnames(X_nons)]
              X_nons <- cbind(y_nons, X_nons[, !grepl("Intercept", colnames(X_nons)), drop=FALSE])

              results_ipw_b_combined <- nonprob_ipw(selection = reformulate(dr_coefs_sel[[1]]),
                                                    target = target,
                                                    data = as.data.frame(X_nons),
                                                    svydesign = svydesign,
                                                    pop_totals = pop_totals_,
                                                    pop_means = pop_means,
                                                    pop_size = pop_size,
                                                    method_selection = method_selection,
                                                    strata = strata,
                                                    case_weights = case_weights[strap_nons],
                                                    na_action = na_action,
                                                    control_selection = control_selection,
                                                    control_inference = control_inference,
                                                    start_selection = start_selection,
                                                    verbose = FALSE,
                                                    se = FALSE,
                                                    pop_size_fixed=pop_size_fixed)

              ## estimate the mi
              results_mi_b_combined <- nonprob_mi(outcome = as.formula(paste0(target_vars, reformulate(dr_coefs_sel[[1]]))),
                                                  data = as.data.frame(X_nons),
                                                  svydesign = svydesign,
                                                  pop_totals = pop_totals_,
                                                  pop_means = pop_means,
                                                  pop_size = pop_size,
                                                  method_outcome = method_outcome,
                                                  family_outcome = family_outcome,
                                                  strata = strata,
                                                  case_weights = case_weights[strap_nons],
                                                  na_action = na_action,
                                                  control_outcome = control_outcome,
                                                  control_inference = control_inference,
                                                  start_outcome = start_outcome,
                                                  verbose = FALSE,
                                                  se = FALSE,
                                                  pop_size_fixed=pop_size_fixed)

              boot_obj[b, ] <- mu_hatDR(y_hat = results_mi_b_combined$output$mean,
                                        y_resid = do.call("cbind", results_mi_b_combined$ys_resid),
                                        weights = case_weights[strap_nons],
                                        weights_nons = results_ipw_b_combined$ipw_weights,
                                        N_nons = sum(results_ipw_b_combined$ipw_weights))

            }
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", b, " iteration: ", e$message, sep = "")
              message(info)
            }
          }
        )
        if (verbose) {
          utils::setTxtProgressBar(pb_boot, b)
        }
        b <- b + 1
      }
    }
  } else {
    # Multicore processing
    if (verbose) message("Multicore bootstrap in progress...")

    cl <- parallel::makeCluster(control_inference$cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl = cl,
                            varlist = NULL,
                            envir = getNamespace("nonprobsvy"))

    if (!is.null(svydesign)) {
      # Parallel bootstrap for probability and non-probability samples
      boot_obj <- foreach::`%dopar%`(
        obj = foreach::foreach(b = 1:num_boot, .combine = rbind),
        ex = {

          strap_rand_svy <- which(rep_weights[, b] != 0)
          weights_rand_strap_svy <- rep_weights[, b] * weights(svydesign)
          pop_size_strap <- sum(weights_rand_strap_svy)

          data_prob <- svydesign$variables[strap_rand_svy, ]
          data_prob$weight <- weights_rand_strap_svy[strap_rand_svy]

          svyd_call <- svydesign$call
          svyd_call$data <- as.name("data_prob")
          svydesign_b <- eval(svyd_call)

          strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / case_weights)

          results_ipw_b <- nonprob_ipw(selection = selection,
                                     target = target,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign_b,
                                     pop_totals = NULL,
                                     pop_means = NULL,
                                     pop_size = NULL,
                                     method_selection = method_selection,
                                     strata = strata,
                                     case_weights = case_weights[strap_nons],
                                     na_action = na_action,
                                     control_selection = control_selection,
                                     control_inference = control_inference,
                                     start_selection = start_selection,
                                     verbose = verbose,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)
          ## estimate the mi
          results_mi_b <- nonprob_mi(outcome = outcome,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign_b,
                                     pop_totals = NULL,
                                     pop_means = NULL,
                                     pop_size = NULL,
                                     method_outcome = method_outcome,
                                     family_outcome = family_outcome,
                                     strata = strata,
                                     case_weights = case_weights[strap_nons],
                                     na_action = na_action,
                                     control_outcome = control_outcome,
                                     control_inference = control_inference,
                                     start_outcome = start_outcome,
                                     verbose = verbose,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)

          if (vars_combine) {

            ipw_coefs_sel <- names(results_ipw_b$selection$coefficients)
            mi_coefs_sel <- lapply(results_mi_b$outcome, coef)
            dr_coefs_sel <- lapply(mi_coefs_sel, function(x) {
              mi_cols <- names(x[abs(x)>0])
              combined <- sort(base::union(ipw_coefs_sel, mi_cols))
              combined[!grepl("Intercept", combined)]
            })

            ## combining variables for selection
            selection_vars <- all.vars(formula.tools::rhs(outcome))
            outcome_vars <- all.vars(formula.tools::rhs(selection))
            target_vars <- all.vars(formula.tools::lhs(outcome))
            combined_vars <- reformulate(union(selection_vars, outcome_vars))

            y_nons <- subset(data[strap_nons, ], select=target_vars)
            X_nons <- model.matrix(combined_vars, data[strap_nons, ])
            X_rand <- model.matrix(combined_vars, svydesign_b$variables) ## if design is present
            X_all <- rbind(X_rand, X_nons)

            X_nons <- cbind(y_nons, X_nons[, !grepl("Intercept", colnames(X_nons)), drop=FALSE])
            X_rand <- X_rand[, !grepl("Intercept", colnames(X_rand)), drop=FALSE]
            svydesign_b_ <- svydesign_b
            svydesign_b_$variables <- cbind(svydesign_b_$variables, X_rand)

            if (bias_corr) {

              print("\nBias correction...")

              ## consider different start
              par0 <- numeric(NCOL(X_all)*2)
              names(par0) <- rep(colnames(X_all), times = 2)

              bias_corr_result_b <- nleqslv::nleqslv(
                x = par0,
                fn = u_theta_beta_dr,
                method = control_selection$nleqslv_method,
                global = control_selection$nleqslv_global,
                xscalm = control_selection$nleqslv_xscalm,
                jacobian = TRUE,
                control = list(
                  scalex = rep(1, length(par0)),
                  maxit = control_selection$maxit
                ),
                R = results_ipw_b$R,
                X = X_all,
                y = c(rep(0, sum(results_ipw_b$R==0)), y_nons[, 1]),
                weights = c(weights(svydesign_b_), results_ipw_b$case_weights),#c(weights(svydesign_), weights),
                method_selection = method_selection,
                family_outcome = family_outcome
              )

              theta_hat <- bias_corr_result_b$x[1:NCOL(X_all)]
              beta_hat <- bias_corr_result_b$x[(NCOL(X_all) + 1):(2 * NCOL(X_all))]

              bias_corr_ps <- method$make_link_inv(unname(drop(X_all %*% theta_hat)))
              bias_corr_ipw_weights <- 1/bias_corr_ps[results_ipw_b$R == 1]
              bias_corr_mu_rand_pred <- as.vector(get(family_outcome)()$linkinv(X_all[results_ipw_b$R == 0, ] %*% beta_hat))
              bias_corr_mu_nons_pred <- as.vector(get(family_outcome)()$linkinv(X_all[results_ipw_b$R == 1, ] %*% beta_hat))
              bias_corr_mu_resid <- bias_corr_mu_nons_pred - y_nons

              boot_obj_b <- mu_hatDR(y_hat = weighted.mean(bias_corr_mu_rand_pred, weights(svydesign_b_)),
                                        y_resid = as.matrix(bias_corr_mu_resid),
                                        weights = results_ipw_b$case_weights,
                                        weights_nons = bias_corr_ipw_weights,
                                        N_nons = sum(bias_corr_ipw_weights))
            } else {
              results_ipw_b_combined <- nonprob_ipw(data = X_nons,
                                                    target = outcome,
                                                    selection =  reformulate(dr_coefs_sel[[1]]),
                                                    svydesign = svydesign_b_,
                                                    pop_totals = NULL,
                                                    pop_means = NULL,
                                                    pop_size = NULL,
                                                    method_selection = method_selection,
                                                    strata = strata,
                                                    case_weights = case_weights[strap_nons],
                                                    na_action = na_action,
                                                    control_selection = control_selection,
                                                    control_inference = control_inference,
                                                    start_selection = start_selection,
                                                    verbose = FALSE,
                                                    se = FALSE,
                                                    pop_size_fixed = pop_size_fixed)
              ## estimate the mi
              results_mi_b_combined <- nonprob_mi(outcome = as.formula(paste0(target_vars, reformulate(dr_coefs_sel[[1]]))),
                                                  data = X_nons,
                                                  svydesign = svydesign_b_,
                                                  pop_totals = NULL,
                                                  pop_means = NULL,
                                                  pop_size = NULL,
                                                  method_outcome = method_outcome,
                                                  family_outcome = family_outcome,
                                                  strata = strata,
                                                  case_weights = case_weights[strap_nons],
                                                  na_action = na_action,
                                                  control_outcome = control_outcome,
                                                  control_inference = control_inference,
                                                  start_outcome = start_outcome,
                                                  verbose = FALSE,
                                                  se = FALSE,
                                                  pop_size_fixed=pop_size_fixed)

              boot_obj_b <- mu_hatDR(y_hat = results_mi_b_combined$output$mean,
                                        y_resid = do.call("cbind", results_mi_b_combined$ys_resid),
                                        weights = case_weights[strap_nons],
                                        weights_nons = results_ipw_b_combined$ipw_weights,
                                        N_nons = sum(results_ipw_b_combined$ipw_weights))

            }
          }


          as.matrix(boot_obj_b)
        }
      )
    } else {
      # Parallel bootstrap for non-probability samples only
      boot_obj <- foreach::`%dopar%`(
        obj = foreach::foreach(b = 1:num_boot, .combine = rbind),
        ex = {

          strap_nons <- sample.int(replace = TRUE, n = NROW(data), prob = 1 / case_weights)

          results_ipw_b <- nonprob_ipw(selection = selection,
                                     target = target,
                                     data = data[strap_nons, ],
                                     svydesign = svydesign,
                                     pop_totals = pop_totals,
                                     pop_means = pop_means,
                                     pop_size = pop_size,
                                     method_selection = method_selection,
                                     strata = strata,
                                     case_weights = case_weights[strap_nons],
                                     na_action = na_action,
                                     control_selection = control_selection,
                                     control_inference = control_inference,
                                     start_selection = start_selection,
                                     verbose = verbose,
                                     se = FALSE,
                                     pop_size_fixed=pop_size_fixed)
          ## estimate the mi
          results_mi_b <- nonprob_mi(outcome = outcome,
                                     data = data[strap_nons, ],
                                   svydesign = svydesign,
                                   pop_totals = pop_totals,
                                   pop_means = pop_means,
                                   pop_size = pop_size,
                                   method_outcome = method_outcome,
                                   family_outcome = family_outcome,
                                   strata = strata,
                                   case_weights = case_weights[strap_nons],
                                   na_action = na_action,
                                   control_outcome = control_outcome,
                                   control_inference = control_inference,
                                   start_outcome = start_outcome,
                                   verbose = verbose,
                                   se = FALSE,
                                   pop_size_fixed=pop_size_fixed)

          ## combination of variables after bootstrap
          if (vars_combine) {

          }
          boot_obj_b <- mu_hatDR(y_hat = results_mi_b$output$mean,
                                 y_resid = do.call("cbind", results_mi_b$ys_resid),
                                 weights = case_weights[strap_nons],
                                 weights_nons = results_ipw_b$ipw_weights,
                                 N_nons = sum(results_ipw_b$ipw_weights))
          as.matrix(boot_obj_b)
        }
      )
    }
  }

  return(boot_obj)
}
