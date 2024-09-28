#
# MATH3871 Lab 2
#

# Question 1

inverse_CDF_1 <- function(y) {
    return((-1 + sqrt(1+8*y))/2)
}

fx_1 <- function(x) {
    return (x + 0.5)
}

obs_1 <- inverse_CDF_1(runif(5000))
hist(obs_1, prob=TRUE)

xx_1 <- seq(0,1,length=100)
print(xx_1)
lines(xx_1, fx_1(xx_1), col="red")

# Question 2

# @precondition - assumes that 0 <= y <= 1
inverse_CDF_2 <- function(y) {
    ifelse(
        y < 0.5,
        0.5 * log(2 * y),
        -0.5 * log(2*(1-y))
    )
}

fx_2 <- function(x) {
    ifelse(
        x < 0,
        exp(2*x),
        exp(-2*x)
    )
}

obs_2 <- inverse_CDF_2(runif(5000))
hist(obs_2, prob=TRUE)
xx_2 <- seq(-6, 6, length=100)
lines(xx_2, fx_2(xx_2), col='red')

# Question 3

inverse_CDF_3 <- function(y, lam) {
    return(-lam * log(1-y))
}

fx_3 <- function(x, lam) {
    return(1/lam * exp(-x/lam))
}

draw_exp_samples <- function(n, lam) {
    return(inverse_CDF_3(runif(n), lam))
}

lam <- 3
xx_3 <- seq(0,25,length=100)
par(mfrow = c(2,2))
for (n in c(50, 250, 500, 1000)) {
    hist(draw_exp_samples(n, lam), prob=TRUE, main=sprintf("Sample size: %d", n))
    lines(xx_3, fx_3(xx_3, lam), col='red')
}
par(mfrow = c(1,1))

# Question 4

# Function to draw a single observation from the Gamma(n, lambda) distribution.
draw_gamma_obs <- function(n, lam) {
    return(sum(draw_exp_samples(n, lam)))
}

# Draw 'num_gamma_samples' samples from Gamma(n, lambda)
lam <- 3
n <- 1000
num_gamma_samples <- 2000
gamma_observations <-numeric(0)
for (i in 1:num_gamma_samples) {
    gamma_observations <- c(gamma_observations, draw_gamma_obs(n, lam))
}

hist(gamma_observations, prob=TRUE)
xx_4 <- seq(2500, 3500, length=100)
lines(xx_4, dgamma(xx_4, shape=n, rate=1/lam), col='red', lwd=3)

