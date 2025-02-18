# Bootstrap for the MI estimator
#' @importFrom survey as.svrepdesign
#' @importFrom nleqslv nleqslv
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar

boot_mi <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    family_outcome,
                    start_outcome,
                    num_boot,
                    weights_rand,
                    mu_hat,
                    svydesign,
                    model_obj = model_obj,
                    rep_type,
                    method,
                    control_outcome,
                    control_inference,
                    pop_totals,
                    verbose,
                    ...) { # TODO add methods instead of conditions

  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  k <- 1
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = num_boot, style = 3)
  }


  pmm_match_type <- control_outcome$pmm_match_type
  nn_exact_se <- control_inference$nn_exact_se
  pmm_reg_engine <- control_outcome$pmm_reg_engine
  pi_ij <- control_inference$pi_ij
  comp2_stat <- numeric(length = num_boot)

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
    if (class(svydesign)[1] != "pps") {
      rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
    } else {
      stop("The pps bootstrap variance estimator is under development")
    }
    if (method == "glm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            # X_rand_strap <- X_rand[which(rep_weights[,k] != 0),]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta <- X_rand %*% beta
            y_strap_rand <- family_nonprobsvy$linkinv(eta)

            # mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
            mu_hat_boot <- weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "nn") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            # strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand)
            # weights_rand_strap <- weights_rand[strap_rand]
            # X_rand_strap <- X_rand[strap_rand, , drop = FALSE]
            # N_strap <- sum(weights_rand_strap)

            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            X_rand_strap <- X_rand[which(rep_weights[, k] != 0), ]
            weights_rand_strap <- weights_rand_strap_svy[strap_rand_svy]

            model_rand <- nonprob_mi_nn(
              data = X_nons_strap,
              query = X_rand_strap,
              k = control_outcome$k,
              treetype = control_outcome$treetype,
              searchtype = control_outcome$searchtype
            )

            y_rand_strap <- apply(model_rand$nn.idx, 1,
              FUN = function(x) mean(y_strap[x])
              # FUN=function(x) mean(sample_nonprob$short_[x])
            )

            mu_hat_boot <- weighted.mean(x = y_rand_strap, w = weights_rand_strap)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "pmm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]
            # using svy package
            strap_rand_svy <- which(rep_weights[, k] != 0)
            weights_rand_strap_svy <- rep_weights[, k] * weights_rand
            N_strap <- sum(weights_rand_strap_svy)
            X_rand_strap <- X_rand[which(rep_weights[, k] != 0), ]
            n_rand_strap <- nrow(X_rand_strap)
            weights_rand_strap <- weights_rand_strap_svy[strap_rand_svy]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta_rand <- X_rand_strap %*% beta
            eta_nons <- X_nons_strap %*% beta
            y_rand_strap <- family_nonprobsvy$linkinv(eta_rand)
            y_nons_strap <- family_nonprobsvy$linkinv(eta_nons)

            model_rand <- switch(control_outcome$pmm_match_type,
              { # 1
                nonprob_mi_nn(
                  data = y_strap,
                  query = y_rand_strap,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              },
              { # 2
                nonprob_mi_nn(
                  data = y_nons_strap,
                  query = y_rand_strap,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              }
            )

            y_rand_strap <- apply(model_rand$nn.idx, 1,
              FUN = function(x) mean(y_strap[x])
              # FUN=function(x) mean(sample_nonprob$short_[x])
            )

            mu_hat_boot <- weighted.mean(x = y_rand_strap, w = weights_rand_strap)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              # info <- paste("iteration ", k, "/", num_boot, ", estimated mean = ", mu_hat_boot, sep = "")
              # print(info)
              utils::setTxtProgressBar(pb, k)
            }
            # slower option
            # if (pmm_exact_se) {
            #   comp2 <- pmm_exact(pi_ij,
            #                      weights_rand,
            #                      n_nons = n_nons,
            #                      y = y,
            #                      pmm_reg_engine = pmm_reg_engine,
            #                      model_obj = model_obj,
            #                      svydesign = svydesign,
            #                      pmm_match_type = pmm_match_type,
            #                      k = control_inference$k,
            #                      N = N)
            #   comp2_stat[k] <- comp2
            # }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta <- pop_totals %*% beta / N
            y_strap_rand <- family_nonprobsvy$linkinv(eta)

            # mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap_svy, N_strap)
            mu_hat_boot <- as.vector(y_strap_rand)
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "nn") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_rand <- nonprob_mi_nn(
              data = X_nons_strap,
              query = t(pop_totals / N),
              k = control_outcome$k,
              treetype = control_outcome$treetype,
              searchtype = control_outcome$searchtype
            )
            mu_hat_boot <- mean(y_strap[model_rand$nn.idx])
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    } else if (method == "pmm") {
      while (k <= num_boot) {
        tryCatch(
          {
            strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
            weights_strap <- weights[strap]
            X_nons_strap <- X_nons[strap, , drop = FALSE]
            y_strap <- y[strap]

            model_strap <- stats::glm.fit(
              x = X_nons_strap,
              y = y_strap,
              weights = weights_strap,
              family = family,
              start = start_outcome
            )

            beta <- model_strap$coefficients
            eta_rand <- pop_totals %*% beta / N
            eta_nons <- X_nons_strap %*% beta
            y_strap_rand <- family_nonprobsvy$linkinv(eta_rand)
            y_strap_nons <- family_nonprobsvy$linkinv(eta_nons)


            model_rand <- switch(control_outcome$pmm_match_type,
              { # 1
                nonprob_mi_nn(
                  data = y_strap,
                  query = y_strap_rand,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              },
              { # 2
                nonprob_mi_nn(
                  data = y_strap_nons,
                  query = y_strap_rand,
                  k = control_outcome$k,
                  treetype = control_outcome$treetype,
                  searchtype = control_outcome$searchtype
                )
              }
            )
            mu_hat_boot <- mean(y_strap[model_rand$nn.idx])
            mu_hats[k] <- mu_hat_boot
            if (verbose) {
              utils::setTxtProgressBar(pb, k)
            }
            k <- k + 1
          },
          error = function(e) {
            if (verbose) {
              info <- paste("An error occurred in ", k, " iteration: ", e$message, sep = "")
              print(info)
            }
          }
        )
      }
    }
  }
  if (method == "pmm") {
    if (nn_exact_se) {
      comp2 <- pmm_exact(pi_ij,
        weights_rand,
        n_nons = n_nons,
        y = y,
        pmm_reg_engine = pmm_reg_engine,
        model_obj = model_obj,
        svydesign = svydesign,
        pmm_match_type = pmm_match_type,
        k = control_inference$k,
        N = N
      )
      comp2 <- mean(comp2_stat)
    } else {
      comp2 <- 0
    }
  } else {
    comp2 <- 0
  }
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2) + comp2
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats,
    comp2 = comp2
  )
}


# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
boot_mi_multicore <- function(X_rand,
                              X_nons,
                              weights,
                              y,
                              family_outcome,
                              start_outcome,
                              num_boot,
                              weights_rand,
                              mu_hat,
                              svydesign,
                              method,
                              control_outcome,
                              control_inference,
                              pop_totals,
                              cores,
                              verbose,
                              ...) {
  # mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  family <- family_outcome
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  rep_type <- control_inference$rep_type

  if (is.character(family_outcome)) {
    family_nonprobsvy <- paste(family_outcome, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  if (verbose) message("Multicores bootstrap in progress...")

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl = cl, varlist = c(
    "model_ps", "start_fit", "est_method_ipw", "control_sel",
    "mle", "mu_hatIPW", "nonprob_mi_nn"
  ))

  if (is.null(pop_totals)) {
    n_rand <- nrow(X_rand)
    N <- sum(weights_rand)
    if (class(svydesign)[1] != "pps") {
      rep_weights <- survey::as.svrepdesign(svydesign, type = rep_type, replicates = num_boot)$repweights$weights
    } else {
      stop("The pps bootstrap variance estimator is under development.")
    }
    if (method == "glm") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta <- X_rand %*% beta
          y_strap_rand <- family_nonprobsvy$linkinv(eta)
          weighted.mean(x = y_strap_rand, w = weights_rand_strap_svy)
        }
      )
    } else if (method == "nn") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[strap_rand_svy]

          model_rand <- nonprob_mi_nn(
            data = X_nons_strap,
            query = X_rand_strap,
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          y_rand_strap <- apply(model_rand$nn.idx, 1,
            FUN = function(x) mean(y_strap[x])
          )
          weighted.mean(x = y_rand_strap, w = weights_strap_rand)
        }
      )
    } else if (method == "pmm") {
      k <- 1:num_boot
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = k, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          # using svy package
          strap_rand_svy <- which(rep_weights[, k] != 0)
          weights_rand_strap_svy <- rep_weights[, k] * weights_rand
          N_strap <- sum(weights_rand_strap_svy)
          X_rand_strap <- X_rand[strap_rand_svy, , drop = FALSE]
          weights_strap_rand <- weights_rand_strap_svy[weights_rand_strap_svy != 0]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )


          beta <- model_strap$coefficients
          eta_rand <- X_rand_strap %*% beta
          eta_nons <- X_nons_strap %*% beta
          y_rand_strap <- family_nonprobsvy$linkinv(eta_rand)
          y_nons_strap <- family_nonprobsvy$linkinv(eta_nons)


          model_rand <- switch(control_outcome$pmm_match_type,
            { # 1
              nonprob_mi_nn(
                data = y_strap,
                query = y_rand_strap,
                k = control_outcome$k,
                treetype = control_outcome$treetype,
                searchtype = control_outcome$searchtype
              )
            },
            { # 2
              nonprob_mi_nn(
                data = y_nons_strap,
                query = y_rand_strap,
                k = control_outcome$k,
                treetype = control_outcome$treetype,
                searchtype = control_outcome$searchtype
              )
            }
          )

          y_rand_strap <- apply(model_rand$nn.idx, 1,
            FUN = function(x) mean(y_strap[x])
          )
          weighted.mean(x = y_rand_strap, w = weights_strap_rand)
        }
      )
    }
  } else {
    N <- pop_totals[1]
    if (method == "glm") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta <- pop_totals %*% beta / N
          y_strap_rand <- family_nonprobsvy$linkinv(eta)
          as.vector(y_strap_rand)
        }
      )
    } else if (method == "nn") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_rand <- nonprob_mi_nn(
            data = X_nons_strap,
            query = t(pop_totals / N),
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          mean(y_strap[model_rand$nn.idx])
        }
      )
    } else if (method == "pmm") {
      mu_hats <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:num_boot, .combine = c),
        ex = {
          strap <- sample.int(replace = TRUE, n = n_nons, prob = 1 / weights)
          weights_strap <- weights[strap]
          X_nons_strap <- X_nons[strap, , drop = FALSE]
          y_strap <- y[strap]

          model_strap <- stats::glm.fit(
            x = X_nons_strap,
            y = y_strap,
            weights = weights_strap,
            family = family,
            start = start_outcome
          )

          beta <- model_strap$coefficients
          eta_rand <- pop_totals %*% beta
          eta_nons <- X_nons_strap %*% beta
          y_strap_rand <- family_nonprobsvy$linkinv(eta_rand)
          y_strap_nons <- family_nonprobsvy$linkinv(eta_nons)


          model_rand <- nonprob_mi_nn(
            data = y_strap_nons,
            query = y_strap_rand,
            k = control_outcome$k,
            treetype = control_outcome$treetype,
            searchtype = control_outcome$searchtype
          )
          mean(y_strap[model_rand$nn.idx])
        }
      )
    }
  }
  # mu_hat_boot <- mean(mu_hats)
  boot_var <- 1 / (num_boot - 1) * sum((mu_hats - mu_hat)^2)
  list(
    var = boot_var,
    # mu = mu_hat_boot,
    stat = mu_hats
  )
}
