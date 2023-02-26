#' controlSel
#'
#' control function for the selection equation in the nonprob function
#' @param method estimation method
#' @param epsilon - a
#' @param maxit - a
#' @param trace - a
#' @param lambda - a
#' @param optim_method - a
#' @param overlap - a
#' @param dependence - a
#'
#' @export

controlSel <- function(method = "glm.fit", #perhaps another control function for model with variables selection
                       epsilon = 1e-8,
                       maxit = 25,
                       trace = FALSE,
                       lambda = 1.25,
                       optim_method = "NR",
                       overlap = FALSE,
                       dependence = FALSE
                       ) {

  list(epsilon = epsilon,
              maxit = maxit,
              trace = trace,
              lambda = lambda,
              optim_method = optim_method,
              overlap = overlap,
              dependence = dependence)

}

#' controlOut
#'
#' control function for the outcome equation in the nonprob function
#' @param method estimation method
#' @param epsilon - a
#' @param maxit - a
#' @param trace - a
#' @param lambda - a
#'
#' @export

controlOut <- function(method = "glm.fit", #perhaps another control function for model with variables selection
                       epsilon = 1e-6,
                       maxit = 25,
                       trace = FALSE,
                       lambda_theta = 1.25,
                       lambda_beta = 0.25
                       ) {

  list(epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       lambda_theta = lambda_theta,
       lambda_beta = lambda_beta)


}


#' controlInf
#'
#' control function for the inference method in the nonprob function
#' @param est_method estimation method
#' @param var_method variance method
#' @param alpha - a
#'
#' @export

controlInf <- function(est_method = c("likelihood",
                                      "integrative"),
                       var_method = c("analytic",
                                      "bootstrap"),
                       type_rep = c("auto", "JK1", "JKn", "BRR", "bootstrap",
                                    "subbootstrap","mrbbootstrap","Fay"),
                       alpha = 0.05) {

  list(est_method = if(missing(est_method)) "likelihood" else est_method,
       var_method = if(missing(var_method)) "analytic" else var_method,
       type_rep = if(missing(type_rep)) "subbootstrap" else type_rep,
       alpha = alpha)

}
