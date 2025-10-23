#################
# MATH3871 Lab 10
#################

# Question 1

# a)
set.seed(3871)
n <- 100
x <- rnorm(n, mean = 0, sd = 1)
y <- 2 + 3*x + rnorm(n, mean = 0, sd = 2)

# Linear regression model
model <- lm(y ~ x)
summary(model)

# Plot the data and regression line
plot(x, y, main = "Linear Regression", xlab = "x", ylab = "y")
abline(model, col = "red")

# b)
# Bayesian linear regression
library(MCMCpack)
bayes_model <- MCMCregress(y ~ x, burnin = 1000, mcmc = 10000)
summary(bayes_model)

# Plot posterior distributions
plot(bayes_model)

# Question 2

# a)
# Simulate data from a mixture model
set.seed(3871)
n <- 200
z <- rbinom(n, 1, 0.7)
x <- rep(NA, n)
x[z == 0] <- rnorm(sum(z == 0), mean = -2, sd = 1)
x[z == 1] <- rnorm(sum(z == 1), mean = 2, sd = 1.5)

# Plot the histogram
hist(x, breaks = 30, freq = FALSE, main = "Mixture Model", xlab = "x")

# b)
# EM algorithm for Gaussian mixture model
em_mixture <- function(x, K = 2, max_iter = 100, tol = 1e-6) {
  n <- length(x)
  
  # Initialize parameters
  pi <- rep(1/K, K)
  mu <- seq(min(x), max(x), length.out = K)
  sigma <- rep(sd(x), K)
  
  # Initialize log-likelihood
  ll_old <- -Inf
  
  for (iter in 1:max_iter) {
    # E-step: Calculate responsibilities
    r <- matrix(NA, n, K)
    for (k in 1:K) {
      r[, k] <- pi[k] * dnorm(x, mu[k], sigma[k])
    }
    r <- r / rowSums(r)
    
    # M-step: Update parameters
    for (k in 1:K) {
      pi[k] <- mean(r[, k])
      mu[k] <- sum(r[, k] * x) / sum(r[, k])
      sigma[k] <- sqrt(sum(r[, k] * (x - mu[k])^2) / sum(r[, k]))
    }
    
    # Calculate log-likelihood
    ll_new <- sum(log(rowSums(sapply(1:K, function(k) pi[k] * dnorm(x, mu[k], sigma[k])))))
    
    # Check convergence
    if (abs(ll_new - ll_old) < tol) break
    ll_old <- ll_new
  }
  
  return(list(pi = pi, mu = mu, sigma = sigma, r = r, iter = iter))
}

# Run EM algorithm
em_result <- em_mixture(x)
print(em_result)

# Plot the fitted mixture model
hist(x, breaks = 30, freq = FALSE, main = "Fitted Mixture Model", xlab = "x")
curve(em_result$pi[1] * dnorm(x, em_result$mu[1], em_result$sigma[1]) + 
      em_result$pi[2] * dnorm(x, em_result$mu[2], em_result$sigma[2]), 
      add = TRUE, col = "red", lwd = 2)
