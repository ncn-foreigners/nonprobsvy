#' @import mathjaxr
NULL
#' @title genSimData
#' @author Łukasz Chrostowski, Maciej Beręsewicz
#
#' @description generate simulated data according to Chen, Li & Wu (2020), section 5
#'
#' \loadmathjax
#
#' @param N \code{integer}, population size, default 10000
#' @param n \code{integer}, big data sample, default 1000
#' @return \code{genSimData} returns a data.frame, with the following columns:
#' \itemize{
#' \item{x0 -- intercept}
#' \item{x1 -- the first variable based on z1}
#' \item{x2 -- the second variable based on z2 and x1}
#' \item{x3 -- the third variable based on z3 and x1 and x2}
#' \item{x4 -- the third variable based on z4 and x1, x2 and x3}
#' \item{y30 -- y generated from the model y=2+x1+x2+x3+x4+sigma*epsilon, so the cor(y,y_hat) = 0.30}
#' \item{y60 -- y generated from the model y=2+x1+x2+x3+x4+sigma*epsilon, so the cor(y,y_hat) = 0.60}
#' \item{y80 -- y generated from the model y=2+x1+x2+x3+x4+sigma*epsilon, so the cor(y,y_hat) = 0.80}
#' \item{rho -- true propensity scores for big data such that sum(rho)=n}
#' \item{srs -- probabilities of inclusion to random sample such that max(srs)/min(srs)=50}
#' }
#'
#' @references Chen, Y., Li, P., & Wu, C. (2020). Doubly Robust Inference With Nonprobability Survey Samples. Journal of the American Statistical Association, 115(532), 2011–2021. https://doi.org/10.1080/01621459.2019.1677241
#'
#' @examples
#' ## generate data with N=20000 and n=2000
#' genSimData(N=20000,n=2000)
#'
#' ## generate data when big data is almost as N
#' genSimData(N=10000,n=9000)
#'
#' @importFrom stats cor lm.fit rbinom rchisq rexp rnorm runif uniroot
#'
#' @export
genSimData <- function(N=10000, n=1000) {

  find_sigma_regsim <- function(sigma, rho) {
    y <- 2 + x1 + x2 + x3 + x4 + sigma*epsilon
    model_fit <- lm.fit(x = cbind(1, x1,x2,x3,x4), y = y)
    res <- cor(y, model_fit$fitted.values) - rho
    res
  }

  find_theta_logsim <- function(theta, n) {
    eta <- theta + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4
    rho <- exp(eta)  / (1 + exp(eta))
    res <- sum(rho) - n
    res
  }

  z1 <- rbinom(N, 1, 0.5)
  z2 <- runif(N,0,2)
  z3 <- rexp(N)
  z4 <- rchisq(N, 4)
  epsilon <- rnorm(N)

  x0 <- rep(1, N)
  x1 <- z1
  x2 <- z2 + 0.3*x1
  x3 <- z3 + 0.2*(x1 + x2)
  x4 <- z4 + 0.1*(x1 + x2 + x3)
  sigma30 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.3)$root
  sigma60 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.6)$root
  sigma80 <- uniroot(find_sigma_regsim, lower = 0, upper = 30, rho = 0.8)$root
  y30 <- 2 + x1 + x2 + x3 + x4 + sigma30*epsilon
  y60 <- 2 + x1 + x2 + x3 + x4 + sigma60*epsilon
  y80 <- 2 + x1 + x2 + x3 + x4 + sigma80*epsilon
  theta <- uniroot(find_theta_logsim, lower = -100, upper = 100, n = n)$root
  rho <-  exp(theta + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4)  / (1 + exp(theta + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4))

  p <- uniroot(f = function(x) max(x3+x) - 50*min(x3+x), lower = -200, upper = 100)$root

  sim_df <- data.frame(x0, x1,x2,x3,x4,y30,y60,y80,rho, srs = (x3+p)/10)

  return(sim_df)

}
