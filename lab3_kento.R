#
# MATH3871 Lab 3
#

### Question 1 ################################################################

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
points(accepted_x, accepted_y, col='blue', pch='.')
points(rejected_x, rejected_y, col='red', pch='.')

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

print_result <- function(simulated, actual) {
    sprintf("Monte Carlo: %f, Actual: %f",
            simulated, actual)
}

# (a)

f <- function(x) {exp(exp(x))}
print_result(mean(f(runif(5000, 0, 1))), integrate(f, 0, 1)$value)

# (b)

f <- function(x) {exp(-x^2)}
u <- runif(5000, 0, 1)
simulated <- 2 * mean(exp(-(1/u - 1)^2) * 1/u^2)
actual <- integrate(f, -Inf, Inf)$value
print_result(simulated, actual)

# (c)

f <- function(x,y) {exp((x+y)^2)}
u1 <- runif(5000,0,1)
u2 <- runif(5000,0,1)
simulated <- mean(exp((u1 + u2)^2))
# dunno how to double integral in R haha
simulated

# (d)

f <- function(x,y) {exp(-(x+y))}
u1 <- runif(5000, 0, 1)
u2 <- runif(5000, 0, 1)
simulated <- mean(exp(-(1/u1 + 1/u2 - 2)) * 1/(u1^2 * u2^2) * (u2 < u1))
simulated

#### Question 3 ################################################################

h <- function(x) {
    return ((cos(50*x) + sin(20*x))^2)
}

sample_vars <- numeric(0)
estimates <- numeric(0)
sample_sizes <- seq(100, 20000, 100)
for (n in sample_sizes) {
    values <- h(runif(n, 0, 1))
    estimate <- mean(values)
    estimates <- c(estimates, estimate)
    sample_var <- 1/n^2 * sum((values - estimate)^2)
    sample_vars <- c(sample_vars, sample_var)
}


# Is this the plot they want?
upper <- estimate + sample_vars
lower <- estimate - sample_vars

plot(c(sample_sizes, sample_sizes), c(upper, lower), cex=0.5, col="lightblue")
abline(h=estimate, lwd=2, col="blue")
text(x=7000, y=0.967, labels = "Monte Carlo estimate (n=200k)", col="blue")

# ...or this?
upper <- estimates + sample_vars
lower <- estimates - sample_vars

plot(
    c(sample_sizes, sample_sizes), c(upper, lower),
    cex=0.5, pch=3, col="lightblue"
)
points(sample_sizes, estimates, cex=0.5, col="blue")
abline(h=0.965, col="red", lwd=2)
text(x=2000, y=0.975, labels = "True value", col="red")


#### Question 4 ################################################################

beta_dist <- function(a, b, x) {
    return (1/beta(a,b) * x^(a-1) * (1-x)^(b-1))
}

xx <- seq(0, 1, length=100)
# Informative prior
plot(xx, beta_dist(27.24, 75.96, xx), type="l", col="#3399FF", lwd=2)
lines(xx, beta_dist(2.24, 0.96, xx), col="#99CCFF", lwd=2)

# Uninformative prior
lines(xx, beta_dist(1, 1, xx), type="l", col="#FFCC33", lwd=2)
lines(xx, beta_dist(26, 76, xx), col="#FFCC00", lwd=2)

# the posterior distributions are quite close... will there be a big difference
# in acceptance rate with rejection sampling?









#### Rejection sampling again

# Target distribution
f <- function(x) {
    return (3/2 * x^3 + 11/8 * x^2 + 1/6 * x + 1/12)
}
# Proposal distribution - Beta(1,1)
g <- function(x) {
    return (dbeta(x, 1, 1))
}

curve(f, 0, 1, xlab='x', ylab='density', lwd=2, col='red') # Target density
xx <- seq(0,1,length=100)
K = max(f(xx)) # K = 3.125
lines(xx, g(xx), col='lightblue', lwd=2) # g(x)
lines(xx, K*g(xx), col='blue', lwd=2) # Kg(x)
text(0.1, 1.2, 'g(x)')
text(0.1, 2.9, 'Kg(x)')

MAX_ITER = 50000
NUM_TO_SAMPLE = 25000

X = rbeta(MAX_ITER, 1, 1) # Sampling values from our proposal, Beta(1,1)
U = runif(MAX_ITER) # This will always be uniform 0,1 - doesn't depend on the
                # proposal distribution
i = 1
accepted = c()
while (i <= MAX_ITER & length(accepted) < NUM_TO_SAMPLE) {
    test_u = U[i]
    test_x = f(X[i]) / (K * dunif(X[i], 0, 1))
    if (test_u <= test_x) {
        # accept
        accepted = c(accepted, X[i])
    }
    i = i + 1
}

hist(accepted, prob=TRUE)
lines(xx, f(xx), lwd=2, col='red')
