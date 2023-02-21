#' Internal functions

bootMI <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   ...
                   ){


  mu_hats <- vector(mode = "numeric", length = num_boot)
  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  N <- sum(weights_rand)

  k <- 1

  while (k <= num_boot) {

    strap_rand <- sample.int(replace = TRUE, n = n_rand, prob = 1/weights_rand) # to change -> replicate weights based on sampling design
                                                                                # see Wu and Thompson (2020) [p. 248]

    h <- c()
    for (x in strap_rand) {

      r <- length(strap_rand[strap_rand == x])
      h <- append(h, r)

    }

    strap <- sample.int(replace = TRUE, n = n_nons)
    weights_strap <- weights[strap]
    Xnons_strap <- X_nons[strap,]
    y_strap <- y[strap]
    weights_rand_strap <- h * weights_rand[strap_rand]
    Xrand_strap <- X_rand[strap_rand, ]
    Nstrap <- sum(weights_rand_strap)


    model_strap <- nonprobMI_fit(x = Xnons_strap,
                                 y = y_strap,
                                 weights = weights_strap,
                                 family_outcome = family_outcome)

    beta <- model_strap$coefficients

    ystrap_nons <- as.numeric(Xnons_strap %*% beta)

    ystrap_rand <- as.numeric(Xrand_strap %*% beta)

    mu_hat_boot <- mu_hatMI(ystrap_rand, weights_rand_strap, N)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mean(mu_hats))^2) # mean may be replaced by mu_hatMI

  boot_var

}

bootIPW <- function(X_rand,
                    X_nons,
                    weights,
                    y,
                    family_outcome,
                    num_boot,
                    weights_rand,
                    mu_hat,
                    dependency,
                    ...){





}

bootDR <- function(X_rand,
                   X_nons,
                   weights,
                   y,
                   family_outcome,
                   num_boot,
                   weights_rand,
                   mu_hat,
                   dependency,
                   ...){



}
