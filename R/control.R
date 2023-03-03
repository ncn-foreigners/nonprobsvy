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
                       epsilon = 1e-6,
                       maxit = 25,
                       trace = FALSE,
                       optim_method = "NR",
                       overlap = FALSE,
                       dependence = FALSE,
                       h_x = c(1, 2)
                       ) {

  list(epsilon = epsilon,
      maxit = maxit,
      trace = trace,
      optim_method = optim_method,
      overlap = overlap,
      dependence = dependence,
      h_x = if(missing(h_x)) 1 else h_x
      )

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

controlOut <- function(method = c("glm", "nn"), #perhaps another control function for model with variables selection
                       epsilon = 1e-6,
                       maxit = 25,
                       trace = FALSE,
                       lambda = 0.25,
                       k = 5
                       ) {

  list(method = if(missing(method)) "glm" else method,
       epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       lambda = lambda,
       k = k)


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
                       bias_inf = c("union", "div"),
                       penalty = c("SCAD", "LASSO"),
                       alpha = 0.05) {

  list(est_method = if(missing(est_method)) "likelihood" else est_method,
       var_method = if(missing(var_method)) "analytic" else var_method,
       type_rep = if(missing(type_rep)) "subbootstrap" else type_rep,
       bias_inf = if(missing(bias_inf)) "union" else bias_inf,
       penalty = if(missing(penalty)) "SCAD" else penalty,
       alpha = alpha)

}
