#' Function for the mass imputation model using pmm method
#'
#' @importFrom stats glm.fit
#' @importFrom ncvreg cv.ncvreg
#' @importFrom stats update
#' @importFrom survey svymean
#'
#' @description
#' Model for the outcome for the mass imputation estimator
#'
#'
#' @param y_nons target variable from non-probability sample
#' @param X_nons a `model.matrix` with auxiliary variables from non-probability sample
#' @param X_rand a `model.matrix` with auxiliary variables from non-probability sample
#' @param weights case / frequency weights from non-probability sample
#' @param svydesign a svydesign object
#' @param family_outcome family for the glm model
#' @param start_outcome start parameters
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals population totals from the `nonprob` function
#' @param pop_size population size from the `nonprob` function
#' @param control_outcome controls passed by the `control_out` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#'
#' @returns an `nonprob_method` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model_fitted}{fitted model either an `glm.fit` or `cv.ncvreg` object}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample}
#'   \item{y_rand_pred}{predicted values for the probability sample or population totals}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{y_mi_hat}{estimated population mean for the target variable}
#'   \item{vars_selection}{whether variable selection was performed}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sampl component}
#'   \item{model}{model type (character `"pmm"`)}
#' }
#' @export
method_pmm <- function(y_nons,
                       X_nons,
                       X_rand,
                       weights,
                       svydesign,
                       family_outcome,
                       start_outcome,
                       vars_selection,
                       pop_totals,
                       pop_size,
                       control_outcome,
                       verbose,
                       se) {
  return(
    structure(
      list(
        model_fitted = NULL,
        y_nons_pred = NULL,
        y_rand_pred = NULL,
        coefficients = NULL,
        svydesign = NULL,
        y_mi_hat = NULL,
        vars_selection = NULL,
        var_prob = NULL,
        var_nonprob = NULL,
        model = "pmm"
      ), class = "nonprob_method")
  )
}
