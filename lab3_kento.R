#
# MATH3871 Lab 3
#

#### Question 1 ################################################################

# This is the distribution to sample from
f <- function(x) {
    return(x^2 * exp(-x))
}

# Normalised distribution
f_norm <- function(x) {
    return(f(x) / integrate(f, lower=0, upper=1)$value)
}

# This is the enveloping function
g <- function(x) {
    return(0.4 * x)
}

xx <- seq(0, 1, length=100)

plot(xx, f(xx), type='l', lwd=2,
     xlab='x', ylab='density', main='Rejection sampling')

lines(c(0,1), c(exp(-1),exp(-1)), lwd=2, col='red')

K <- exp(-1)
N <- 5000
sim_x <- runif(N, 0, 1)
sim_y <- runif(N, 0, K)

accepted_x <- sim_x[sim_y < f(sim_x)]
accepted_y <- sim_y[sim_y < f(sim_x)]
rejected_x <- sim_x[sim_y >= f(sim_x)]
rejected_y <- sim_y[sim_y >= f(sim_x)]

# See accepted vs rejected points
points(accepted_x, accepted_y, col='blue')
points(rejected_x, rejected_y, col='red')

# Check histogram of simulated values vs distribution function
hist(accepted_x, prob=TRUE)
lines(xx, dgamma(xx, 3, 1)/pgamma(1,3,1), lwd=3, col='blue')

# Notes
# - This isn't the best way to do it because it doesn't ensure N samples from
#   the distribution. Should make a function that repeats individual simulations
#   until N samples are accepted.
# - If g(x) is used to sample values for x*, then Kg(x) must be used to sample
#   values for y* - i.e. g(x) has to stay the same.


#### Question 2 ################################################################

# (a)

sprintf("Monte Carlo: %f", mean(exp(exp(runif(5000, 0, 1)))))
f <- function(x) {exp(exp(x))}
sprintf("Actual: %f", integrate(f, 0, 1)$value)

# (b)






