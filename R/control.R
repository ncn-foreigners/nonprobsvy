#' controlSel
#'
#' control function for the selection equation in the nonprob function
#' @param method estimation method
#' @param epsilon - a
#' @param maxit - a
#' @param trace - a
#' @param lambda - a
#'
#' @export

controlSel <- function(method = "glm.fit", #perhaps another control function for model with variables selection
                       epsilon = 1e-8,
                       maxit = 25,
                       trace = FALSE,
                       lambda = 1.25
                       ) {

  return(list(epsilon = epsilon,
              maxit = maxit,
              trace = trace,
              lambda = lambda
  ))

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
                       epsilon = 1e-8,
                       maxit = 25,
                       trace = FALSE,
                       lambda = 0.25
                       ) {

  return(list(epsilon = epsilon,
              maxit = maxit,
              trace = trace,
              lambda = lambda
              ))

}


#' controlInf
#'
#' control function for the inference method in the nonprob function
#' @param est.method estimation method
#' @param var.method variance method
#' @export

controlInf <- function(est.method = c("likelihood",
                                      "integrative"),
                       var.method = c("analytic",
                                      "bootstrap")) {

  list(est.method = if(missing(est.method)) "likelihood" else est.method,
       var.method = if(missing(var.method)) "analytic" else est.method)

}
