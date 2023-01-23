#' controlSel
#'
#' control function for the selection equation in the nonprob function
#' @param method estimation method
#'
#' @export
controlSel <- function(method = "glm.fit"
                       #epsilon = 1e-8,
                       #maxit = 25,
                       #trace = FALSE
                       ) {

}

#' controlOut
#'
#' control function for the outcome equation in the nonprob function
#' @param method estimation method
#' @export
controlOut <- function(method = "glm.fit",
                       epsilon = 1e-8,
                       maxit = 25,
                       trace = FALSE
                       ) {

  return(list(epsilon = 1e-8,
              maxit = 25,
              trace = FALSE))

}


#' controlInf
#'
#' control function for the inference method in the nonprob function
#' @param est.method estimation method
#' @param var.method variance method
#' @export

controlInf <- function(est.method, var.method) {

}
