cv_nonprobsvy <- function(X, R, weights_X, method_selection, h, nfolds = 10) {

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  X_nons <- cbind(X[loc_nons,], weights_X[loc_nons], R[loc_nons])
  X_rand <- cbind(X[loc_rand,], weights_X[loc_rand], R[loc_nons])
  k <- 1

  loss_theta_av <- vector(mode = "numeric", length = 50)
  lambdas <- setup_lambda(X = X,
                          y = R,
                          weights = weights_X,
                          method_selection = method_selection,
                          lambda_min = 0,
                          nlambda = 50)


  X_nons <- X_nons[sample(nrow(X_nons)), ]
  X_rand <- X_rand[sample(nrow(X_rand)), ]

  folds_nons <- cut(seq(1,nrow(X_nons)), breaks=nfolds, labels=FALSE) #split nonprobabability sample into K parts
  folds_rand <- cut(seq(1,nrow(X_rand)), breaks=nfolds, labels=FALSE) #split probabability sample into K parts

  # pair K subsets randomly
  sample_nons <- sample(1:nfolds, nfolds, replace = FALSE)
  sample_rand <- sample(1:nfolds, nfolds, replace = FALSE)

  for (lambda in lambdas) {

    loss_theta_vec <- vector(mode = "numeric", length = nfolds)

    for(i in 1:nfolds){

      # train data for X_nons
      idx_nons <- which(folds_nons==sample_nons[i], arr.ind=TRUE)
      X_nons_train <- X_nons[-idx_nons, ]


      # test data for X_nons
      X_nons_test <- X_nons[idx_nons, ]


      # train data for X_rand
      idx_rand <- which(folds_rand==sample_rand[i], arr.ind=TRUE)
      X_rand_train <- X_rand[-idx_rand, ]

      # test data for X_rand
      X_rand_test <- X_rand[idx_rand, ]


      X_train <- rbind(X_rand_train, X_nons_train)
      X_test <- rbind(X_rand_test, X_nons_test)


      theta_est <- do.call("fit", X_train, X_train$R,  X_test$weights_X, method_selection, h)

      loss_theta_vec[i] <- loss_theta(par = theta_est, R = X_test$R, X = X_test,
                                      weights = X_test$weights_X, h = h, method_selection = method_selection)


    }

    loss_theta_av[k] <- mean(loss_theta_vec)

    k <- k + 1

  }

  i <- which.min(loss_theta_av)

  print(lambdas[i])

  lambdas[i]

}



setup_lambda <- function(X, y, weights, method_selection, alpha = 1, lambda_min, log,lambda = FALSE, nlambda, ...) { #consider panalty factor here

  fit <- glm.fit(x = X, y = y, weights = weights, family = binomial(link = method_selection))
  #fit <- glm(y~1, weights = weights, family = binomial(link = method_selection))

  n <- length(y)
  p <- ncol(X)
  w <- fit$weights
  r <- residuals(fit, "working") * w
  zmax <- .Call("maxprod", X, r)/n
  lambda.max <- zmax/alpha


  if (log.lambda) { # lambda sequence on log-scale
    if (lambda.min==0) {
      lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
    } else {
      lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
    }
  } else { # lambda sequence on linear-scale
    if (lambda.min==0) {
      lambda <- c(seq(lambda.max, 0.001*lambda.max, length = nlambda-1), 0)
    } else {
      lambda <- seq(lambda.max, lambda.min*lambda.max, length = nlambda)
    }
  }
  lambda
}


fit <- function(X, R, weights, method_selection, h) {

  p <- ncol(X)
  init_theta <- rep(0, p)
  # variables selection using score equation for theta
  par0 <- init_theta
  LAMBDA <- Matrix::Matrix(matrix(0, p, p), sparse = TRUE)
  it <- 0
  for(jj in 1:maxit) {
    it <- it + 1
    if (it == maxit) break

    u_theta0 <- u_theta(par = par0, R = R, X = X,
                        weights = weights_X, h = h,
                        method_selection = method_selection)

    u_theta0_der <- u_theta_der(par = par0, R = R, X = X,
                                weights = weights_X, h = h,
                                method_selection = method_selection)

    diag(LAMBDA) <- abs(q_lambda(par0, lambda))/(eps + abs(par0))
    par <- par0 + MASS::ginv(as.matrix(u_theta0_der + LAMBDA)) %*% (u_theta0 - LAMBDA %*% par0) # perhaps 'solve' function instead of 'ginv'
    # equation (13) in article
    if (sum(abs(par - par0)) < eps) break;
    if (sum(abs(par - par0)) > 1000) break;

    par0 <- par

  }

  par <- as.vector(par)
  theta_est <- par
  theta_est

}
