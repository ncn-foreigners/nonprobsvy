#' nonprobOv
#'
#' nonprobOv: Function for correcting selection bias considering the possible overlapping and dependency between the nonprobability and probability sample, function needs furher work
#'
#' @importFrom stats glm.fit
#' @importFrom betareg betareg.fit


# consider function structure as in nonprobDR, nonprobMI, nonprobIPW, nonprobSel

nonprobOv <- function(X_nons,
                      X_rand,
                      d_rand,
                      dependent,
                      method_selection) {

  betamodel <- betareg.fit(x = X_rand,
                           y = 1/d_rand,
                           link = method_selection)

  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  inv_link <- method$make_link_inv
  prob <- inv_link(X_nons %*% betamodel$coefficients$mean)
  d_rnons <- 1/prob

  # drnons is an inclusion probability in the probability sample for units in the nonprobability sample
  d <- c(d_rnons, d_rand)

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X <- rbind(X_nons, X_rand)
  XY <- cbind(X, R)

  O <- O_hat_model(x = X,
                   X = X,
                   R = R)
  L <- NULL

  if (dependent) {

    L <- L_hat_model(x = X_rand,
                     X = X,
                     R = R_rand)

    ps <- O * L/(d - 1)
    weights <- 1/ps[loc_nons]

  } else {
    ps <- O/(O + d - 1)
    weights <- 1/ps[loc_nons]
  }


  structure(
    list(weights = weights,
         O_hat = O[loc_nons],
         L_hat = L,
         d_rnons = d_rnons)
  )
}


# consider merging overlap and bootstrap functions in one, function needs further work

boot_overlap <- function(X_rand,
                         X_nons,
                         y,
                         weights, #for the nonprobability sample
                         O_hat,
                         L_hat,
                         d_rand,
                         d_rnons,
                         num_boot = 1000,
                         dependency,
                         N) {


  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  mu_hats <- vector(mode = "numeric", length = num_boot)
  d <- c(d_rand, d_rnons)
  R_nons <- rep(1, n_nons)
  R_rand <- rep(0, n_rand)
  R <- c(R_nons, R_rand)
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)


  weights <- weights*N/sum(weights)
  prob_floor <- weights - floor(weights)
  weights_strap <- weights
  weights_idx <- sample.int(length(weights), replace = TRUE, prob = prob_floor)
  weights_strap[weights_idx] <- ceiling(weights_strap[weights_idx])
  weights_strap[!weights_idx] <- floor(weights_strap[!weights_idx])

  pseudo_pop <- X_nons[rep(1:n_nons, times = weights_strap), ] # pseudo population - to consider
  pseudo_d_rnons <- rep(d_rnons, times = weights_strap)
  pseudo_weights <- rep(weights, times = weights_strap)


  k <- 1

  if (!dependency) {


    while (k <= num_boot) {

      strap_rand <- sample.int(n = nrow(pseudo_pop), size = n_rand, replace = TRUE, prob = 1/pseudo_d_rnons)
      strap_nons <- sample.int(n = nrow(pseudo_pop), size = n_nons, replace = TRUE, prob = 1/pseudo_weights)
      #strap_d <- sample.int(n = n_nons + n_rand, replace = TRUE)
      d_rnons_strap <- pseudo_d_rnons[c(strap_rand, strap_nons)]

      X_rand_strap <- pseudo_pop[strap_rand, ]
      X_nons_strap <- pseudo_pop[strap_nons, ]
      X_strap <- rbind(X_nons_strap, X_rand_strap)
      R_nons_strap <- rep(1, nrow(X_nons_strap))
      R_rand_strap <- rep(0, nrow(X_rand_strap))
      R_strap <- c(R_nons_strap, R_rand_strap)
      loc_nons_strap <- which(R_strap == 1)
      loc_rand_strap <- which(R_strap == 0)

      # Removing overlapping units here - to do

      O_strap <- O_hat_model(x = X_strap,
                             X = X_strap,
                             R = R_strap)

      ps_strap <- O_strap/(O_strap + d_rnons_strap - 1)
      weights_strap <- 1/ps_strap[loc_nons_strap]
      N_strap <- sum(weights_strap)

      mu_hat_boot <- mu_hatIPW(y = y,
                               weights = weights_strap,
                               N = N_strap)

      mu_hats[k] <- mu_hat_boot

      k <- k + 1

    }

  } else {

    pseudo_L <- rep(L_hat[loc_nons], times = weights_strap)
    pseudo_O <- rep(O_hat[loc_nons], times = weights_strap)

    # L <- L_hat_model(x = pseudo_pop,
                 #    X = X,
                  #   R = pseudo_R)

    ps1 <- 1 - pseudo_L
    ps0 <- 1 - pseudo_d_rnons + pseudo_d_rnons*pseudo_O*pseudo_L + pseudo_L*(pseudo_d_rnons - 1)


    # weights for strap_nons to compute

    while (k <= num_boot) {

      strap_rand <- sample.int(n = nrow(pseudo_pop), size = n_rand, replace = TRUE, prob = 1/pseudo_d_rnons)

      ps1[!strap_rand] <- ps0
      strap_nons <- sample.int(n = nrow(pseudo_pop), size = n_nons, replace = TRUE, prob = 1/ps1)
      d_rnons_strap <- pseudo_d_rnons[c(strap_rand, strap_nons)]

      X_rand_strap <- pseudo_pop[strap_rand, ]
      X_nons_strap <- pseudo_pop[strap_nons, ]
      X_strap <- rbind(X_nons_strap, X_rand_strap)
      R_nons_strap <- rep(1, nrow(X_nons_strap))
      R_rand_strap <- rep(0, nrow(X_rand_strap))
      R_strap <- c(R_nons_strap, R_rand_strap)
      loc_nons_strap <- which(R_strap == 1)
      loc_rand_strap <- which(R_strap == 0)


      O_strap <- O_hat_model(x = X_strap,
                             X = X_strap,
                             R = R_strap)


      L_strap <- L_hat_model(x = X_rand_strap,
                             X = X_strap,
                             R = R_rand_strap)

      ps_strap <- O_strap * L_strap/(d_rnons_strap - 1)
      weights_strap <- 1/ps_strap[loc_nons_strap]
      N_strap <- sum(weights_strap)

      mu_hat_boot <- mu_hatIPW(y = y,
                               weights = weights_strap,
                               N = N_strap)

      mu_hats[k] <- mu_hat_boot

      k <- k + 1

    }

  }

  var <- 1/(num_boot-1) * sum((mu_hats - mean(mu_hats))^2)

  var

}


L_hat_model <- function(x,
                        X,
                        R) {

  modelL <- stats::glm.fit(x = x,
                           y = R,
                           family = binomial(link = "logit")) # glm.fit: algorithm did not converge
                                                              # to consider

  L <- 1/(1 + exp(X %*% modelL$coefficients))

  L

}



O_hat_model <- function(x,
                        X,
                        R) {

  modelO <- stats::glm.fit(x = x,
                           y = R,
                           family = binomial(link = "logit"))


  O <- exp(X %*% modelO$coefficients)

  O

}


