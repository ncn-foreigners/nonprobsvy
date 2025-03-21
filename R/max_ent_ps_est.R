#' An internal function to compute maximum entropy propensity score estimator
#'
#' @description
#' Fits a maximum entropy propensity score estimator for non-probability surveys
#' based on a density ratio model. The estimator is obtained by solving the
#' moment condition:
#' \deqn{\sum_{i\in A}\delta_i\exp(\phi_0+\phi_1'x_i)(x_i-\hat{\bar{x}}_0)=0,}
#' where \eqn{\hat{\bar{x}}_0} is the weighted covariate mean for the probability sample,
#' with \eqn{\hat{N}_0=\sum_{i\in A}d_i-\sum_{i\in A}\delta_i}.
#'
#' @param X Matrix of covariates from the probability sample.
#' @param delta Binary indicator vector; 1 indicates membership in the non-probability sample.
#' @param d Numeric vector of design weights for the probability sample.
#' @param start_phi1 Numeric vector of starting values for \eqn{\phi_1} (default is a zero vector).
#' @param maxit Integer, maximum number of iterations for the optimization (default 100).
#' @param tol Convergence tolerance for the iterative solver (default 1e-6).
#' @param nleqslv_method Character, specifies the method for \code{nleqslv} optimization (default "Broyden").
#' @param nleqslv_global Character, global strategy for \code{nleqslv} (default "dbldog").
#' @param nleqslv_xscalm Character, scaling method for \code{nleqslv} (default "auto").
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{phi}: A list with components \code{phi0} and \code{phi1} (the estimated parameters).
#'   \item \code{theta_PS2}: The final maximum entropy propensity score estimator.
#'   \item \code{diagnostics}: The output from the \code{nleqslv} solver.
#' }
#'
#' @details
#' The estimator is based on the model:
#' \deqn{\log\{r(x)\}=\phi_0+\phi_1'x,}
#' which, after a reparameterization, is equivalent to
#' \deqn{\pi(x)=\frac{1}{1+\exp(\phi^*_0+\phi_1'x)},}
#' with \eqn{\phi^*_0=\log(\pi_0/\pi_1)+\phi_0}.
#'
#' Using Horvitz-Thompson estimation, the sample-based moment condition is
#' formulated as
#' \deqn{\sum_{i\in A}\delta_i\exp(\phi_0+\phi_1'x_i)(x_i-\hat{\bar{x}}_0)=0,}
#' where
#' \deqn{\hat{\bar{x}}_0=\frac{1}{\hat{N}_0}\left\{\sum_{i\in A}d_ix_i-\sum_{i\in A}\delta_ix_i\right\}, \quad \hat{N}_0=\sum_{i\in A}d_i-\sum_{i\in A}\delta_i.}
#'
#' The intercept \eqn{\phi_0} is determined by enforcing
#' \deqn{\sum_{i\in A}\delta_i\exp(\phi_0+\phi_1'x_i)=\sum_{i\in A}\delta_i,}
#' which leads to
#' \deqn{\phi_0=\log\left\{\sum_{i\in A}\delta_i\right\}-\log\left\{\sum_{i\in A}\delta_i\exp(\phi_1'x_i)\right\}.}
#'
#' Finally, the propensity score estimator is computed as
#' \deqn{\hat{\theta}_{PS2}=\frac{1}{\sum_{i\in A}d_i}\sum_{i\in A}\delta_i\left\{1+\frac{\hat{N}_0}{\sum_{i\in A}\delta_i}\exp(\phi_0+\phi_1'x_i)\right\}.}
#'
#' @importFrom nleqslv nleqslv
#'
#' @noRd
est_method_maxent <- function(X, delta, d,
                              start_phi1 = rep(0, ncol(X)),
                              maxit = 100,
                              tol = 1e-6,
                              nleqslv_method = "Broyden",
                              nleqslv_global = "dbldog",
                              nleqslv_xscalm = "auto") {
  
  
  X <- as.matrix(X)
  p <- ncol(X)
  N_total <- sum(d)
  N1 <- sum(delta)
  N0_hat <- sum(d) - N1
  
  # hat{x0} = (sum d_i x_i - sum delta_i x_i) / (sum d_i - sum delta_i)
  x_bar0 <- (colSums(d * X) - colSums(X * delta)) / N0_hat
  
  # f(phi1) = sum_{i in A} delta_i * exp(phi0 + phi1' x_i) * (x_i - x_bar0)
  # with phi0 computed as: phi0 = log(N1) - log(sum_{i in A} delta_i * exp(phi1' x_i))
  moment_fun <- function(phi1) {
    eta <- as.vector(X %*% phi1)
    sum_exp <- sum(delta * exp(eta))
    phi0 <- log(N1) - log(sum_exp)
    moments <- colSums((delta * exp(phi0 + eta)) * (X - matrix(rep(x_bar0, each = nrow(X)),
                                                               nrow = nrow(X),
                                                               byrow = FALSE)))
    return(moments)
  }
  
  # Solve for phi1 using nleqslv
  sol <- nleqslv::nleqslv(x = start_phi1, fn = moment_fun,
                          method = nleqslv_method,
                          global = nleqslv_global,
                          xscalm = nleqslv_xscalm,
                          control = list(maxit = maxit, ftol = tol))
  
  phi1_hat <- sol$x
  eta_hat <- as.vector(X %*% phi1_hat)
  sum_exp_hat <- sum(delta * exp(eta_hat))
  phi0_hat <- log(N1) - log(sum_exp_hat)
  
  # Compute the final maximum entropy propensity score estimator:
  propensity_component <- 1 + (N0_hat / N1) * exp(phi0_hat + eta_hat)
  theta_PS2 <- sum(delta * propensity_component) / N_total
  
  list(
    phi = list(phi0 = phi0_hat, phi1 = phi1_hat),
    theta_PS2 = theta_PS2,
    diagnostics = sol
  )
}


# Example 1: Book example
set.seed(123)
N <- 10000          
p <- 10              

X_pop <- matrix(rnorm(N * p, mean = 2, sd = 1), ncol = p)
colnames(X_pop) <- paste0("x", 1:p)


delta_pop <- ifelse(X_pop[, 1] > 2, 1, 0)

prob_idx <- which(delta_pop == 0)
n_A <- 500
sample_A_idx <- sample(prob_idx, n_A)


N_prob <- length(prob_idx)
w_A <- rep(N_prob / n_A, n_A)

sample_B_idx <- which(delta_pop == 1)
w_B <- rep(0, length(sample_B_idx))


indices <- c(sample_A_idx, sample_B_idx)
X_sample <- X_pop[indices, ]
delta_sample <- delta_pop[indices].
d_sample <- numeric(length(indices))
d_sample[1:length(sample_A_idx)] <- w_A
d_sample[(length(sample_A_idx) + 1):length(indices)] <- 0

res_ex1 <- est_method_maxent(X_sample, delta_sample, d_sample)
print(res_ex1)


# Example 2: Paper example  
set.seed(456)
N <- 10000  

region    <- sample(1:17, N, replace = TRUE)
age_group <- sample(1:14, N, replace = TRUE)
sex       <- sample(0:1, N, replace = TRUE)
height    <- rnorm(N, mean = 165, sd = 10)
weight    <- rnorm(N, mean = 70, sd = 15)
waist     <- rnorm(N, mean = 80, sd = 10)
alcohol   <- rbinom(N, 1, 0.3)
systolic  <- rnorm(N, mean = 120, sd = 15)
diastolic <- rnorm(N, mean = 80, sd = 10)
fasting   <- rnorm(N, mean = 90, sd = 10)
creatinine<- rnorm(N, mean = 1, sd = 0.2)


X_pop2 <- cbind(region, age_group, sex, height, weight, waist, alcohol, systolic, diastolic, fasting, creatinine)
colnames(X_pop2) <- c("region", "age_group", "sex", "height", "weight", "waist",
                      "alcohol", "systolic", "diastolic", "fasting", "creatinine")

logit_p <- -3 + 0.1 * height - 0.05 * weight + 0.2 * sex
prob_inclusion <- exp(logit_p) / (1 + exp(logit_p))
delta_pop2 <- rbinom(N, 1, prob_inclusion)


n_A2 <- 1092
sample_A_idx2 <- sample(1:N, n_A2)

w_A2 <- rep(N / n_A2, n_A2)


sample_B_idx2 <- which(delta_pop2 == 1)
w_B2 <- rep(0, length(sample_B_idx2)) 

indices2 <- c(sample_A_idx2, sample_B_idx2)
X_sample2 <- X_pop2[indices2, ]
delta_sample2 <- delta_pop2[indices2]
d_sample2 <- numeric(length(indices2))
d_sample2[1:length(sample_A_idx2)] <- w_A2
d_sample2[(length(sample_A_idx2) + 1):length(indices2)] <- 0


res_ex2 <- est_method_maxent(X_sample2, delta_sample2, d_sample2)
print(res_ex2)

