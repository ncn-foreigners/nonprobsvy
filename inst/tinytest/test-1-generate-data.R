library(sampling)
library(survey)

N <- 10000
n_A <- 500
p <- 50
alpha_vec1 <- c(-2, 1, 1, 1,1, rep(0, p-5))
alpha_vec2 <- c(0,0,0,3,3,3,3, rep(0, p-7))
beta_vec <- c(1,0,0,1,1,1,1, rep(0, p-7))

## generate X
X <- cbind("(Intercept)"=1, matrix(rnorm(N*(p-1)), nrow=N, byrow=T,
                                   dimnames = list(NULL, paste0("X",1:(p-1)))))
X_formula  <- as.formula(paste("~", paste0("X",1:(p-1), collapse = " + ")))
X_totals <- colSums(X)
X_means <- colMeans(X[,-1])

## generate Y
Y_11 <- 1 + as.numeric(X %*% beta_vec) +   rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)

## generate probs
pi_A <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A) ## inclusion based on Y_11 only
pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability
pi_B2 <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) -
                  sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## PSM II: nonlinear

## generate
flag_B1 <- rbinom(N, 1, prob = pi_B1)
flag_B2 <- rbinom(N, 1, prob = pi_B2)
flag_A <- UPpoisson(pik = pi_A)

## generate data
pop_data <- data.frame(pi_A, flag_A, flag_B1, flag_B2, Y_11, Y_12, Y_21, Y_22, X[, 2:p])

## generate nonprob data
sample_B1 <- subset(pop_data, flag_B1 == 1)
sample_B2 <- subset(pop_data, flag_B2 == 1)

## generate prob data
sample_A <- subset(pop_data, flag_A == 1)

## population totals
X_totals <- colSums(X)
X_means <- colMeans(X[,-1])

## probability sample
sample_A_svy <- svydesign(ids = ~ 1,
                          probs = ~ pi_A,
                          pps = poisson_sampling(sample_A$pi_A),
                          data = sample_A)
sample_A_svy_cal <- calibrate(sample_A_svy,
                              formula = as.formula(paste0("~", paste(names(X_totals)[2:p], collapse = "+"))),
                              population = X_totals,
                              calfun = cal.raking)
