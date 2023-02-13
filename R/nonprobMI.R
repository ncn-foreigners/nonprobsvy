#' nonprobMI
#
#' nonprobMI: Function for inference based on nonprobability surveys using mass imputation
#
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param method_outcome - a `character` with method for response variable estimation
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @export

nonprobMI <- function(outcome,
                      data,
                      svydesign,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      weights,
                      na_action,
                      control_outcome = controlOut(),
                      control_inference = controlInf(),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  y_name <- colnames(XY_nons)[1]
  X_nons <- model.matrix(XY_nons, data) # matrix of the  nonprobability sample
  svydesign$variables[,y_name] <- rep(0, nrow(svydesign$variables))
  X_rand <- model.matrix(outcome, svydesign$variables) # matrix of the probability sample
  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  N_est_rand <- sum(d_rand)


  ## estimation
  model_nons <- nonprobMI_fit(x = X_nons,
                              y = y_nons,
                              weights = weights,
                              family_outcome = family_outcome)

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  # updating probability sample by adding y_hat variable
  svydesign <- stats::update(svydesign,
                             .y_hat_MI = y_rand_pred)


  nonprobMI_inference <- function(...) {

    mu_hat <- mu_hatMI(y = y_rand_pred,
                       weights = d_rand,
                       N = N_est_rand)


    n_nons <- nrow(X_nons)
    n_rand <- nrow(X_rand)


    # design based variance estimation based on approximations of the second-order inclusion probabilities

    s <- y_rand_pred
    ci <- n_rand/(n_rand-1) * (1 - ps_rand)
    B_hat <- sum(ci*(s/ps_rand))/sum(ci)
    ei <- (s/ps_rand) - B_hat
    db_var <- sum(ci*(ei^2))

    v_b <- 1/N_est_rand^2 * db_var  # first component

    # v_b <- 1/N_estB * t(model$coefficients) * E * model$coefficients

    mx <- 1/N_est_rand * colSums(d_rand * X_rand)

    mh <- 0
    for (i in 1:n_nons) { # matrix product instead of a loop in a near future

      xx <- t(X_nons[i,]) %*% X_nons[i,]
      mh <- mh + xx
    }

    c <- 1/(1/n_nons * mh) %*% mx
    e <- XY_nons[, 1] - y_nons_pred

    # second component
    v_a <- 1/n_nons^2 * t(as.matrix(e^2))  %*% (as.matrix(X_nons) %*% t(as.matrix(c)))^2

    # variance
    var <- v_a + v_b

    se <- sqrt(var)

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    ci <- c(mu_hat - z*se, mu_hat + z*se)

    # bootstrap variance
    boot_var <- bootMI(X_rand,
                       X_nons,
                       weights,
                       y_nons,
                       family_outcome,
                       1000,
                       d_rand,
                       n_nons,
                       n_rand,
                       mu_hat)



    structure(
      list(population_mean = mu_hat,
           variance = var,
           standard_error = se,
           CI = ci,
           beta = model_nons_coefs,
           boot_variance = boot_var
           ),
      class = "Mass imputation")

  }

  # inference based on mi method
  infer_nons <- nonprobMI_inference()

  infer_nons
}


#' nonprobMI_fit
#
#' nonprobMI_fit: Function for outcome variable estimation based on nonprobability sample and using model based approach
#'
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param control_outcome - a
#' @param start - a
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param verbose - a
#' @param model - a
#' @param x - a
#' @param y - a
#'


nonprobMI_fit <- function(outcome,
                          data,
                          weights,
                          svydesign,
                          family_outcome,
                          control_outcome = controlOut(),
                          start = NULL,
                          verbose,
                          model,
                          x,
                          y) {


  family <- family_outcome

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  model_nons <- stats::glm.fit(x = x,
                               y = y,
                               weights = weights,
                               start = start,
                               control = list(control_outcome$epsilon,
                                              control_outcome$maxit,
                                              control_outcome$trace),
                               family = family)


  model_nons

}

#' nonprobMI_nn
#
#' nonprobMI_nn: Function for outcome variable estimation based on nonprobability sample and using predictive mean matching
#'
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param query - a
#' @param k - a
#' @param treetype - a
#' @param searchtype - a
#' @param radius - a
#' @param eps - a


nonprobMI_nn <- function(data,
                         query,
                         k,
                         treetype,
                         searchtype,
                         radius,
                         eps) {


  model_nn <- nn2(data = data,
                  quety = query,
                  k = k,
                  treetype = treetype,
                  searchtype = searchtype,
                  radius = radius,
                  eps = eps)

  model_nn

}

#' mu_hatMI
#
#' mu_hatMI: Function for outcome variable estimation based on mass imputation
#' @param y - a
#' @param weights - a
#' @param N - a

mu_hatMI <- function(y, weights, N) {

  mu_hat <- 1/N * sum(weights * y)

  mu_hat


}


