#
# MATH3871 Lab 4
#

### Question 1 ################################################################

# (a)

# When the question says "draw samples...using importance sampling", it doesn't
# actually mean draw iid samples from the target distribution - it just wants us
# to get samples from the proposal distribution along with their weights, and to
# calculate the ESS.

# Let proposal be uniform 0,1 (as in Lab 3)
n = 5000
samples = runif(n)

# Compute weights by f(x_i) / g(x_i)
weights = samples^2 * exp(-samples) / 1
weights_norm = weights / sum(weights)

ess_unif = 1/(sum(weights_norm^2))
sprintf("ESS from %d samples is %.2f", n, ess_unif)

# (b)

# Redoing rejection sampling for practice...
n = 5000
K = exp(-1)

accepted = c()
i = 1
while (length(accepted) < n) {
    x = runif(1)
    if (runif(1) <= (x^2*exp(-x)/K)) { # with probability f(x*) / Kg(x*)..
        accepted = c(accepted, x)
    }
    i = i + 1
}

# We can also resample from our weighted samples to get unweighted samples which
# are (approximately) representative of the target distribution.
resampled = sample(samples, n, replace=T, prob=weights_norm)
resampled

par(mfrow=c(2,1))
hist(accepted, prob=T, breaks=20, main="Rejection sampling")
hist(resampled, prob=T, breaks=20, main="Importance sampling (with resampling)")
par(mfrow=c(1,1))

sprintf("Mean from rejection sampling: %.4f", mean(accepted))
sprintf("Mean from importance sampling: %.4f", sum(samples * weights_norm))
sprintf("Mean from importance sampling (with resampling): %.4f", mean(resampled))

# (c)
n = 5000
a = 2
b = 1
samples = rbeta(n, a, b)
weights = (samples^2 * exp(-samples)) / dbeta(samples, a, b)
weights_norm = weights / sum(weights)
ess_beta = 1/sum(weights_norm^2)

resampled = sample(samples, n, replace=T, prob=weights_norm)
hist(resampled, prob=T, main="Importance sampling (with resampling) with Beta proposal")
lines(xx, dbeta(xx, a, b), col="blue", lwd=2) # proposal
lines(xx,dgamma(xx,3,1)/pgamma(1,3,1),col='gray', lwd=2) # target

sprintf("Using Unif(0,1) as proposal distribution gave ESS of %.2f", ess_unif)
sprintf("Using Beta(%d,%d) as proposal distribution gave ESS of %.2f", a, b, ess_beta)
# Higher ESS = more efficient choice of proposal, thus Beta proposal is better
# (loosely, ESS represents the number of samples with equal weight - higher ESS
# means lower variability in weights)

#### Question 2 ################################################################

transformed_f <- function(u) {
    exp(-(1/u - 1)^2) * 1/u^2
}

# Plot transformed_f(u) - remember that the area we want to calculate is 2x the
# area of this
xx <- seq(0,1,length=100)
curve(transformed_f, 0, 1, lwd=2)

# Set proposal dist here
proposal <- function(x) {
    return(x/x) # Uniform
}
lines(xx, proposal(xx), lwd=2, col="red")

n <- 10000
samples <- runif(n)
weights <- transformed_f(samples) / proposal(samples)
weights_norm <- weights/sum(weights) # We only use the normalised weights for
# the ESS, not for the approximation itself - this is because we're not dealing
# with a PDF, just trying to approximate integral

area <- mean(weights) * 2
ess <- 1 / sum(weights_norm^2)

sprintf("Estimated area as %.4f with %d samples and ESS %.2f", area, n, ess)

# Now, approximate integral of original expression without transformation

f <- function(x) {
    return (exp(-x^2))
}
curve(f, -10, 10, lwd=2)
xx <- seq(-10, 10, length=100)

lines(xx, dnorm(xx, 0, 1), lwd=2, col='red')

samples <- rnorm(n, 0, 1)
weights <- f(samples) / dnorm(samples, 0, 1)
weights_norm <- weights / sum(weights)

area <- mean(weights)
area

ess <- 1/sum(weights_norm^2)
sprintf("Estimated area as %.4f with %d samples and ESS %.2f", area, n, ess)

#### Question 3 ################################################################

# (a)

samples <- runif(1000)
weights <- dnorm(samples, 0.5, 0.01) / 1
weights_norm <- weights / sum(weights)
ess <- 1 / sum(weights_norm^2)
ess # ESS ~45 from 1000 samples... this is very inefficient

# (b)

while (ess < 5000) {
    samples <- c(samples, runif(1000))
    weights <- dnorm(samples, 0.5, 0.01) /1
    weights_norm <- weights / sum(weights)
    ess <- 1/sum(weights_norm^2)
}
sprintf("%d samples required to get ESS > 5000 (reached ESS = %.2f)", length(samples), ess)

# (c)

samples <- runif(1000)
weights <- (dnorm(samples, 0.5, 0.01))^0.02 / 1
weights_norm <- weights/sum(weights)

resampled <- sample(samples, 1000, replace=TRUE, prob=weights_norm)

hist(resampled, prob=T)
xx <- seq(0.2, 0.8, length=100)
lines(xx, dnorm(xx, 0.5, 0.01), lwd=3, col="lightblue")
lines(xx, (dnorm(xx, 0.5, 0.01))^0.02, lwd=3, col="blue")
lines(xx, dunif(xx, 0, 1), lwd=3, col="red")

mean <- mean(resampled)
sd <- sd(resampled)

# Now, with the N(mean, sd) as proposal distribution...

ess <- 0
samples <- rnorm(1000, mean, sd)
while (ess < 5000) {
    samples <- c(samples, rnorm(1000, mean, sd))
    weights <- dnorm(samples, 0.5, 0.01) / dnorm(samples, mean, sd)
    weights_norm <- weights / sum(weights)
    ess <- 1/sum(weights_norm^2)
}
sprintf("%d samples required to get ESS > 5000 (reached ESS = %.2f)", length(samples) + 1000, ess)
# Sequential importance sampling is much more efficient!
# Why?
# - breaks down inefficient importance sampling stages into a sequence of more
#   efficient stages in order to control the variability of importance weights
