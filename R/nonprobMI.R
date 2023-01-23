#' nonprobMI
#
#' nonprobMI: Function for inference based on nonprobability surveys using mass imputation
#
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param method.outcome - a `character` with method for response variable estimation
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action a
#' @param control.outcome a
#' @param control.inference a
#' @param start a
#' @param verbose a
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @export

nonprobMI <- function(outcome,
                      data,
                      svydesign,
                      method.outcome,
                      family.outcome,
                      subset,
                      weights,
                      na.action,
                      control.outcome = controlOut(),
                      control.inference = controlInf(),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  XY_nons <- model.frame(outcome, data)
  y_name <- colnames(XY_nons)[1]
  X_nons <- model.matrix(XY_nons, data)
  svydesign$variables[,y_name] <- rep(0, nrow(svydesign$variables))
  X_rand <- model.matrix(outcome, svydesign$variables)
  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  N_est_rand <- sum(d_rand)

  ## estimation
  model_nons <- nonprobMI.fit(x = X_nons,
                              y = y_nons)

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  svydesign <- update(svydesign,
                      .y_hat_MI = y_rand_pred) # updating probability sample by adding y_hat variable


  nonprobMI.inference <- function(...) {


    mu_hat <- 1/N_est_rand * sum(d_rand * y_rand_pred)


    n_nons <- nrow(X_nons)
    n_rand <- nrow(X_rand)


    s <- y_rand_pred # is multiplied by psB?
    ci <- n_rand/(n_rand-1) * (1 - ps_rand)
    B_hat <- sum(ci*(s/ps_rand))/sum(ci)
    ei <- (s/ps_rand) - B_hat
    db_var <- sum(ci*(ei^2))

    v_b <- 1/N_est_rand^2 * db_var  # rechange N_estA on N_estB

    # v_b <- 1/N_estB * t(model$coefficients) * E * model$coefficients

    mx <- 1/N_est_rand * colSums(d_rand * X_rand) # rechange N_estA on N_estB

    suma = 0
    for(i in 1:n_nons){

      c <- t(X_nons[i,]) %*% X_nons[i,]
      suma <- suma + c
    }

    c <- 1/(1/n_nons * suma) %*% mx
    e <- XY_nons[, 1] - y_nons_pred

    v_a <- 1/n_nons^2 * t(as.matrix(e^2))  %*% (as.matrix(X_nons) %*% t(as.matrix(c)))^2

    var <- v_a + v_b

    ci <- c(mu_hat - 1.96*sqrt(var), mu_hat + 1.96*sqrt(var))



    return(list("Population mean estimator" = mu_hat,
                "variance" = var,
                "CI" = ci
                ))

  }

  ## inference based on mi method
  infer_nons <- nonprobMI.inference()

  return(infer_nons)
}


nonprobMI.fit <- function(outcome,
                          data,
                          weights,
                          svydesign,
                          family.outcome,
                          control.outcome,
                          verbose,
                          model,
                          x,
                          y) {


  model_nons <- stats::glm.fit(x = x,
                               y = y,
                               # weights = weights,
                               # start = start,
                               # control = control.outcome,
                               # family = family.outcome
                               )

  return(model_nons)

}


