BootMI <- function(Xrand,
                   Xnons,
                   weights,
                   y,
                   num_boot,
                   d,
                   n.nons,
                   n.rand,
                   mu_hat,
                   ...
                   ){


  mu_hats <- vector(mode = "numeric", length = num_boot)

  k <- 1

  while (k <= num_boot) {

    strap_rand <- sample.int(replace = FALSE, n = n.rand) # to change -> replicate weights based on sampling design
    strap <- sample.int(replace = TRUE, n = n.nons)
    weights.strap <- weights[strap]
    Xnons_strap <- Xnons[strap,]
    y_strap <- y[strap]
    d.strap <- d[strap_rand]
    Nstrap <- sum(d.strap)


    model_strap <- nonprobMI.fit(x = Xnons_strap,
                                 y = y_strap,
                                 weights = weights.strap)

    beta <- model_strap$coefficients

    ystrap_nons <- as.numeric(Xnons_strap %*% beta)

    ystrap_rand <- as.numeric(Xrand %*% beta)

    mu_hat_boot <- mu_hatMI(ystrap_rand, d.strap, Nstrap)
    mu_hats[k] <- mu_hat_boot

    k <- k + 1

  }

  boot_var <- 1/num_boot * sum((mu_hats - mu_hat)^2)

  return(boot_var)

}
