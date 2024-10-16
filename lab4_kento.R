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
