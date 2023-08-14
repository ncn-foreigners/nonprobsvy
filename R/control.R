#' @title Control parameters for selection model
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#' @description \code{controlSel} constructs a list with all necessary control parameters
#' for selection model.
#'
#' \loadmathjax
#'
#' @param method estimation method.
#' @param epsilon Tolerance for fitting algorithms by default \code{1e-6}.
#' @param maxit Maximum number of iterations.
#' @param trace logical value. If `TRUE` trace steps of the fitting algorithms. Default is `FALSE`
#' @param optimizer -
#' @param optim_method maximisation method that will be passed to [stats::optim()] function. Default is `BFGS`.
#' @param maxLik_method maximisation method that will be passed to [maxLik::maxLik()] function. Default is `NR`.
#' @param dependence logical value - `TRUE` if samples are dependent.
#' @param key binary key variable
#' @param est_method_sel -
#' @param ov_method -
#' @param h_x Smooth function for the estimating equations.
#' @param penalty -
#' @param a_SCAD -
#' @param a_MCP -
#' @param lambda A user-specified lambda value.
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .001.
#' @param nlambda The number of lambda values. Default is 50.
#' @param nfolds The number of folds for cross validation. Default is 10.
#' @param print_level -
#'
#' @export

controlSel <- function(method = "glm.fit", #perhaps another control function for model with variables selection
                       epsilon = 1e-6,
                       maxit = 500,
                       trace = FALSE,
                       optimizer = c("maxLik", "optim"),
                       maxLik_method = "NR",
                       optim_method = "BFGS",
                       dependence = FALSE,
                       key = NULL,
                       est_method_sel = c("mle", "gee", "mm"),
                       ov_method = c("ev", "mle", "gee"),
                       h_x = c("1", "2"),
                       penalty = c("SCAD", "lasso", "MCP"),
                       a_SCAD = 3.7,
                       a_MCP = 3,
                       lambda = -1,
                       lambda_min = .001,
                       nlambda = 50,
                       nfolds = 10,
                       print_level = 0
                       ) {

  list(epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       optimizer = if(missing(optimizer)) "optim" else optimizer,
       maxLik_method = maxLik_method,
       optim_method = optim_method,
       dependence = dependence,
       key = key,
       est_method_sel = if(missing(est_method_sel)) "mle" else est_method_sel,
       ov_method = if(missing(ov_method)) "ev" else ov_method,
       h_x = if(missing(h_x)) "1" else h_x,
       penalty = if(missing(penalty)) "SCAD" else penalty,
       a_SCAD = a_SCAD,
       a_MCP = a_MCP,
       lambda_min = lambda_min,
       nlambda = nlambda,
       nfolds = nfolds,
       lambda = lambda,
       print_level = print_level
      )

}

#' @title Control parameters for outcome model
#' @description \code{controlOUT} constructs a list with all necessary control parameters
#' for outcome model.
#' @param epsilon Tolerance for fitting algorithms. Default is \code{1e-6}.
#' @param maxit Maximum number of iterations.
#' @param trace logical value. If `TRUE` trace steps of the fitting algorithms. Default is `FALSE`.
#' @param k The k parameter in the [RANN2::nn()] function. Default is 5.
#' @param penalty penalty algorithm for variable selection. Default is `SCAD`
#' @param a_SCAD -
#' @param a_MCP -
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .001.
#' @param nlambda The number of lambda values. Default is 100.
#' @param nfolds -
#' @param treetype -
#' @param searchtype -
#'
#' @export

controlOut <- function(epsilon = 1e-6,
                       maxit = 100,
                       trace = FALSE,
                       k = 5,
                       penalty = c("SCAD", "lasso", "MCP"),
                       a_SCAD = 3.7,
                       a_MCP = 3,
                       lambda_min = .001,
                       nlambda = 100,
                       nfolds = 10,
                       treetype = "kd",
                       searchtype = "standard"
                       ) {

  list(epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       k = k,
       penalty = if(missing(penalty)) "SCAD" else penalty,
       a_SCAD = a_SCAD,
       a_MCP = a_MCP,
       lambda_min = lambda_min,
       nlambda = nlambda,
       nfolds = nfolds,
       treetype = treetype,
       searchtype = searchtype)

}


#' @title Control parameters for inference
#' @description \code{controlINF} constructs a list with all necessary control parameters
#' for statistical inference.
#' @param vars_selection - .
#' @param var_method variance method.
#' @param rep_type replication type for weights in the bootstrap method for variance estimation. Default is `subbootstrap`.
#' @param bias_inf inference method in the bias minimization. Default is `union`.
#' @param num_boot -
#' @param alpha Significance level, Default is 0.05.
#'
#' @export

controlInf <- function(vars_selection = FALSE,
                       var_method = c("analytic",
                                      "bootstrap"),
                       rep_type = c("auto", "JK1", "JKn", "BRR", "bootstrap",
                                    "subbootstrap","mrbbootstrap","Fay"),
                       bias_inf = c("union", "div"),
                       num_boot = 500,
                       alpha = 0.05) {

  list(vars_selection = if(missing(vars_selection)) FALSE else vars_selection,
       var_method = if(missing(var_method)) "analytic" else var_method,
       rep_type = if(missing(rep_type)) "subbootstrap" else rep_type,
       bias_inf = if(missing(bias_inf)) "union" else bias_inf,
       num_boot = num_boot,
       alpha = alpha)

}
