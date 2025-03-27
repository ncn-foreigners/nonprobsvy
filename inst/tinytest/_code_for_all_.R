## testing bootstrap
set.seed(2024)

# Load required data
data(admin)
data(jvs)

# Create objects ----------------------------------------------------------

# Create survey design object
jvs_svy <- svydesign(
    ids = ~1,
    weights = ~weight,
    strata = ~size + nace + region,
    data = jvs
)

N <- sum(weights(jvs_svy))
pop_totals <- colSums(model.matrix(~region + private + nace + size, jvs)*jvs$weight)
pop_means <- pop_totals[-1]/N


# simulated data (Kim & Yang 2019) ----------------------------------------------------------

kim2019_N <- 1e5 ## 1000000
n <- 500
x1 <- rnorm(n = kim2019_N, mean = 1, sd = 1)
x2 <- rexp(n = kim2019_N, rate = 1)
epsilon <- rnorm(n = kim2019_N)
y1 <- 1 + x1 + x2 + epsilon
y2 <- 0.5*(x1 - 0.5)^2 + x2 + epsilon
p1 <- exp(x2)/(1+exp(x2))
p2 <- exp(-0.5+0.5*(x2-2)^2)/(1+exp(-0.5+0.5*(x2-2)^2))
flag_bd1 <- rbinom(n = kim2019_N, size = 1, prob = p1)
flag_srs <- as.numeric(1:kim2019_N %in% sample(1:kim2019_N, size = n))
base_w_srs <- kim2019_N/n
population <- data.frame(x1,x2,y1,y2,p1,p2,base_w_srs, flag_bd1, flag_srs)
base_w_bd <- kim2019_N/sum(population$flag_bd1)

kim2019_sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs, data = subset(population, flag_srs == 1))
kim2019_sample_nonprob <- subset(population, flag_bd1 == 1)
kim2019_y_true <- c(mean(y1), mean(y2))
kim2019_totals <- colSums(model.matrix(~ x1 + x2, population))


# simulated high-dim data (Yang 2020) -------------------------------------



