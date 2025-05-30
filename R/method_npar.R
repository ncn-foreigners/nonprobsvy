#' Mass imputation using non-parametric model method
#'
#' @importFrom stats update
#' @importFrom survey svymean
#' @importFrom stats loess.control
#' @importFrom stats loess
#' @importFrom stats coef
#'
#' @description
#' Model for the outcome for the mass imputation estimator using loess via `stats::loess`.
#' Estimation of the mean is done using the \eqn{S_B} probability sample.
#'
#' @details Analytical variance
#'
#' The variance of the mean is estimated based on the following approach
#'
#' (a) non-probability part  (\eqn{S_A} with size \eqn{n_A}; denoted as `var_nonprob` in the result)
#'
#' \deqn{
#' \hat{V}_1 = \frac{1}{N^2} \sum_{i=1}^{n_A} \left\lbrace\hat{g}_B(\boldsymbol{x}_i)\right\rbrace^{2} \hat{e}_i^2,
#' }
#'
#' where \eqn{\hat{e}_i=y_i - \hat{m}(x_i)} is the residual and \eqn{\hat{g}_B(\boldsymbol{x}_i) = \left\lbrace \pi_B(\boldsymbol{x}_i) \right\rbrace^{-1}} can be estimated
#' various ways. In the package we estimate \eqn{\hat{g}_B(\boldsymbol{x}_i)} using \eqn{\pi_B(\boldsymbol{x}_i)=E(R | \boldsymbol{x})} as suggested by Chen et al. (2022, p. 6). In particular,
#' we currently support this using stats::loess` with `"gaussian"` family.
#'
#' (b) probability part (\eqn{S_B} with size \eqn{n_B}; denoted as `var_prob` in the result)
#'
#'  This part uses functionalities of the `{survey}` package and the variance is estimated using the following
#'  equation:
#'
#' \deqn{
#' \hat{V}_2=\frac{1}{N^2} \sum_{i=1}^{n_B} \sum_{j=1}^{n_B} \frac{\pi_{i j}-\pi_i \pi_j}{\pi_{i j}}
#' \frac{\hat{m}(x_i)}{\pi_i} \frac{\hat{m}(x_j)}{\pi_j}.
#' }
#'
#' Note that \eqn{\hat{V}_2} in principle can be estimated in various ways depending on the type of the design and whether population size is known or not.
#'
#'
#' @param y_nons target variable from non-probability sample
#' @param X_nons a `model.matrix` with auxiliary variables from non-probability sample
#' @param X_rand a `model.matrix` with auxiliary variables from non-probability sample
#' @param svydesign a svydesign object
#' @param weights case / frequency weights from non-probability sample (default NULL)
#' @param family_outcome family for the glm model)
#' @param start_outcome a place holder (not used in `method_npar`)
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals a place holder (not used in `method_npar`)
#' @param pop_size population size from the `nonprob` function
#' @param control_outcome controls passed by the `control_out` function
#' @param control_inference controls passed by the `control_inf` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#'
#' @returns an `nonprob_method` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model_fitted}{fitted model object returned by `stats::loess`}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample}
#'   \item{y_rand_pred}{predicted values for the probability sample or population totals}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{y_mi_hat}{estimated population mean for the target variable}
#'   \item{vars_selection}{whether variable selection was performed}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sampl component}
#'   \item{model}{model type (character `"npar"`)}
#' }
#'
#' @references
#' Chen, S., Yang, S., & Kim, J. K. (2022). Nonparametric mass imputation for data integration.
#' Journal of Survey Statistics and Methodology, 10(1), 1-24.
#'
#' @examples
#'
#' set.seed(123123123)
#' N <- 10000
#' n_a <- 500
#' n_b <- 1000
#' n_b1 <- 0.7*n_b
#' n_b2 <- 0.3*n_b
#' x1 <- rnorm(N, 2, 1)
#' x2 <- rnorm(N, 2, 1)
#' y1 <- rnorm(N, 0.3 + 2*x1+ 2*x2, 1)
#' y2 <- rnorm(N, 0.3 + 0.5*x1^2+ 0.5*x2^2, 1)
#' strata <- x1 <= 2
#' pop <- data.frame(x1, x2, y1, y2, strata)
#' sample_a <- pop[sample(1:N, n_a),]
#' sample_a$w_a <- N/n_a
#' sample_a_svy <- svydesign(ids=~1, weights=~w_a, data=sample_a)
#' pop1 <- subset(pop, strata == TRUE)
#' pop2 <- subset(pop, strata == FALSE)
#' sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
#'                   pop2[sample(1:nrow(pop2), n_b2), ])
#' res_y_npar <- nonprob(outcome = y1 + y2 ~ x1 + x2,
#'                       data = sample_b,
#'                       svydesign = sample_a_svy,
#'                       method_outcome = "npar")
#' res_y_npar
#'
#' @export
method_npar <- function(y_nons,
                        X_nons,
                        X_rand,
                        svydesign,
                        weights=NULL,
                        family_outcome="gaussian",
                        start_outcome=NULL,
                        vars_selection=FALSE,
                        pop_totals=NULL,
                        pop_size=NULL,
                        control_outcome=control_out(),
                        control_inference=control_inf(),
                        verbose=FALSE,
                        se=TRUE) {

  ## this does not handle aggregated data

  if (missing(y_nons) | missing(X_nons)) {
    stop("`y_nons` and `X_nons`, `X_rand` are required.")
  }

  if (missing(svydesign) | is.null(svydesign)) {
    stop("The NPAR method is suited only for the unit-level data.")
  }

  if (vars_selection) {
    warning("Variable selection for `method_npar` is not yet implemented.
             Estimation is based on the whole set of auxiliary variables.")
    vars_selection <- FALSE
  }

  if (is.null(weights)) weights <- rep(1, nrow(X_nons))
  if (is.null(pop_size)) pop_size <- sum(weights(svydesign))

  if (family_outcome != "gaussian") {
    message("The `method_npar` currently supports only `gaussian` family. Overrwitting to `gaussian`.")
    family_outcome <- "gaussian"
  }

  X <- rbind(X_rand, X_nons)
  R <- rep(c(0, 1), times = c(nrow(X_rand), nrow(X_nons)))

  # NOTE: This is left for the case if we would like to switch to np instead of loess
  # if (verbose) message("Estimation based on `np::npregbw` started")
  # model_fitted <- np::npregbw(xdat = X_nons[,-1], ydat = y_nons, ckertype="gaussian")

  model_fitted <- stats::loess(
    formula = y_nons ~ X_nons[, -1],
    family = family_outcome,
    control = control_outcome$npar_loess)

  y_nons_pred <- predict(model_fitted, X_nons[, -1])
  y_rand_pred <- predict(model_fitted, X_rand[, -1])
  residuals <- as.vector(y_nons - y_nons_pred)
  svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
  svydesign_mean <- survey::svymean( ~ y_hat_MI, svydesign_updated)
  y_mi_hat <- as.numeric(svydesign_mean)


  ## variance components
  if (se) {

      var_prob <- as.vector(attr(svydesign_mean, "var"))

      g_hat <- stats::loess(
        formula = R ~ X[, -1],
        family = family_outcome,
        control = control_outcome$npar_loess)

      var_nonprob <- 1/pop_size^2*sum( g_hat$fitted[R==1]^{-2}*residuals^2)

      var_total <- var_prob + var_nonprob

      # NOTE: This is for further development if we would like to use np instead of loess
      # if (verbose) message("Estimation of variance started. This may take a while...")
      # g_hat_bw <- npcdensbw(ydat = factor(R), xdat = X[,-1], tol = 0.1, ftol = 0.1, nmulti = 2)
      # g_hat <- np::npconmode(bws = g_hat_bw)
      # var_nonprob <- 1/pop_size^2*sum( stats::fitted(g_hat)[R==1]^{-2}*residuals^2)

  } else {
    var_prob <- var_nonprob <- var_total <- NA
  }

  return(
    structure(
      list(
        model_fitted = model_fitted,
        y_nons_pred = y_nons_pred,
        y_rand_pred = y_rand_pred,
        coefficients = coef(model_fitted),
        svydesign = if (is.null(svydesign)) svydesign else svydesign_updated,
        y_mi_hat = y_mi_hat,
        vars_selection = vars_selection,
        var_prob = var_prob,
        var_nonprob = var_nonprob,
        var_total = var_total,
        model = "npar",
        family = family_outcome
      ),
      class = "nonprob_method")
  )
}

