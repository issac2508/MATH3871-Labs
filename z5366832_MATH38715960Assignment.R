#Template Script for Assignment, MATH3871/5960

library(mvtnorm)
library(MASS)
set.seed(1234)

#-------------------------------------------------------------------------------
#INSTRUCTIONS
#1. Do not delete or edit existing text in this template, only add to it.
#2. Ensure your student number is in the file name.

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

########################### PART I Frequentist ##################################
#-------------------------------------------------------------------------------
#QUESTION 1.

wine = read.csv(
    '/Users/kentoseki/Library/CloudStorage/OneDrive-UNSW/UNSW Courses/24T3/MATH3871/Assignment/winequality-red.csv',
    head = TRUE
)

sprintf("Dataset contains %d NA values", sum(is.na(wine)))
# there are 0 NAs

#-------------------------------------------------------------------------------
#QUESTION 2.

wine$good = as.integer(wine$quality >= 6.5)

#-------------------------------------------------------------------------------
#QUESTION 3.

res = glm(
    good ~ fixed.acidity + volatile.acidity + citric.acid + residual.sugar
            + chlorides + free.sulfur.dioxide + total.sulfur.dioxide + density
            + pH + sulphates + alcohol,
    family=binomial,
    data=wine
)

mleest = coef(res)

#-------------------------------------------------------------------------------

############################# PART II Bayesian #################################
#-------------------------------------------------------------------------------
#QUESTION 2.

lpost.LR <- function(beta,X,y) {
    k = ncol(X)

    alpha = rep(0, k)
    omega = 100 * diag(k)
    log.prior = log(dmvnorm(x=beta, mean=alpha, sigma=omega))

    log.likelihood = t(y) %*% X %*% beta - sum(log(1 + exp(X %*% beta)))

    return (log.prior + log.likelihood)
}

#### Just making sure it works
# y = as.matrix(wine$good)
# X = wine[,!names(wine) %in% c("quality", "good")]
# X = as.matrix(cbind("intercept"=rep(1, nrow(X)), X))
# beta = mleest
# lpost.LR(beta, X, y)


#-------------------------------------------------------------------------------
#QUESTION 3.

mhmcmc <- function(y, X, B, nsims, Sigma) {
    # B: the initial beta vector to use

    beta_star = c()  # the proposed value at each iteration
    beta_mat = c()   # the value of the chain at each iteration
    accprob = c()    # the acceptance probability at each iteration
    acc = 0          # the total number of accepted samples

    beta_star[[1]] = B # initialise
    beta_mat[[1]] = B
    accprob[[1]] = NA
    acc = 1

    for (i in 2:nsims) {
        B_prop = mvrnorm(n=1, mu=B, Sigma=Sigma)

        # compute acceptance threshold
        if (exp(lpost.LR(B, X, y)) == 0) {
            log.alpha = 0
        } else {
            log.alpha = min(0, lpost.LR(B_prop, X, y) - lpost.LR(B, X, y))
        }
        log.u = log(runif(1))

        # if accept, update B and append to lists
        if (log.u <= log.alpha) {
            beta_mat[[i]] = B_prop   # accept proposed value
            B = B_prop
            acc = acc + 1
        } else {
            beta_mat[[i]] = B        # reject; keep old value
        }
        beta_star[[i]] = B_prop
        accprob[[i]] = exp(log.alpha)
    }

    # return variables
    mhout = list("beta_mat"=beta_mat, "accprob"=accprob, "beta_star" = beta_star, "acc"=acc)
    return(mhout)
}


#-------------------------------------------------------------------------------
#QUESTION 4.

#Covariance for proposal:
Sigma = diag(c(1100, 0.0015, 0.04, 0.05, 0.0005, 0.4, 0.00001, 0.000001, 0.0001, 0.09, 0.03, 0.002)) #Do not modify

y = as.matrix(wine$good)
X = wine[,!names(wine) %in% c("quality", "good")]
X = as.matrix(cbind("intercept"=rep(1, nrow(X)), X))
B = mleest

nSims = 10^5

mhout1 = mhmcmc(y, X, B, nSims, Sigma)

#(a) produce trace plots:

par(mfrow=c(6,2))
for (i in 1:ncol(X)) {
    # for each ith coefficient...
    param_values = lapply(mhout1$beta_mat, function(l) l[[i]]) # get list of accepted values for the current coefficient
    param_values = unlist(param_values)

    plot(param_values, type="l", main=names(mhout1$beta_mat[[1]])[i], xlab="iteration")
}

# Just doing for a single i to see how to access all the values
# i = 1
# param_values = lapply(mhout1$beta_mat, function(l) l[[i]])
# param_values = unlist(param_values)
#
# par(mfrow=c(1,1))
# plot(param_values, type="l")
# param_values

#(b) proportion of accepted moves:
mh4acc = mhout1$acc / nSims
mh4acc

#-------------------------------------------------------------------------------
#QUESTION 5.

prior.mean = rep(0, ncol(X))
too.high = rep(100, ncol(X))
too.low = rep(-100, ncol(X))

nSims = 10^5

mhout1 = mhmcmc(y, X, mleest, nSims, Sigma)
mhout2 = mhmcmc(y, X, prior.mean, nSims, Sigma)
mhout3 = mhmcmc(y, X, too.high, nSims, Sigma)
mhout4 = mhmcmc(y, X, too.low, nSims, Sigma)

par(mfrow=c(6,2), mar = c(2,2,2,2))
for (i in 1:ncol(X)) {
    # for each ith coefficient...

    param_values_1 = lapply(mhout1$beta_mat, function(l) l[[i]]) # get list of accepted values for the current coefficient
    param_values_1 = unlist(param_values_1)

    param_values_2 = lapply(mhout2$beta_mat, function(l) l[[i]])
    param_values_2 = unlist(param_values_2)

    param_values_3 = lapply(mhout3$beta_mat, function(l) l[[i]])
    param_values_3 = unlist(param_values_3)

    param_values_4 = lapply(mhout4$beta_mat, function(l) l[[i]])
    param_values_4 = unlist(param_values_4)

    ylim = range(param_values_1, param_values_2, param_values_3, param_values_4)

    plot(param_values_1, type="l", main=names(mhout1$beta_mat[[1]])[i], ylim=ylim)
    lines(param_values_2, col="blue")
    lines(param_values_3, col="green")
    lines(param_values_4, col="purple")

}

# individual
# i=2
# param1 = unlist(lapply(mhout1$beta_mat, function(l) l[[i]]))
# param2 = unlist(lapply(mhout2$beta_mat, function(l) l[[i]]))
# param3 = unlist(lapply(mhout3$beta_mat, function(l) l[[i]]))
# param4 = unlist(lapply(mhout4$beta_mat, function(l) l[[i]]))
# ylim = range(param1, param2, param3, param4)
#
# plot(param1, type="l", main=names(mhout1$beta_mat[[1]])[i], ylim=ylim)
# lines(param2, col="blue")
# lines(param3, col="green")
# lines(param4, col="purple")


#-------------------------------------------------------------------------------
#QUESTION 6.

prop = 0.055 # set covariance proportional to MLE covariance
newSigma = vcov(res) * prop
nSims = 10^5

mhout6 = mhmcmc(y, X, B, nSims, newSigma)

mh6acc = mhout6$acc / nSims
mh6acc

#(a)produce trace plots:

par(mfrow=c(6,2), mar=c(2,2,2,2))
for (i in 1:ncol(X)) {
    # for each ith coefficient...
    param_values = lapply(mhout6$beta_mat, function(l) l[[i]])
    param_values = unlist(param_values)

    plot(param_values, type="l", main=names(mhout6$beta_mat[[1]])[i], xlab="iteration")
}

#(b) produce marginal histgrams and overlay MLE:

# Burn in (look at trace plots)
# Remove first 10K samples - all variables look definitely converged by then
converged.samples = tail(mhout6$beta_mat, -10000)

# Check autocorrelation - do we need to do thinning?
par(mfrow=c(6,2))
for (i in 1:ncol(X)) {
    param_values = lapply(converged.samples, function(l) l[[i]])
    param_values = unlist(param_values)
    acf(param_values, lag.max=100)
    # acf(param_values, lag.max=500) # clearly the lag is about 500 samples...
                                     # if we 10^5 / 500 we would get 200 supposedly
                                     # independent samples
}
# Try thinning by taking every 500th sample
# thinned.samples = converged.samples[1:180 * 500]
# Try only taking every 100th sample
# thinned.samples = converged.samples[1:900 * 100]

# Plot histograms
samples = converged.samples # or thinned.samples

par(mfrow=c(6,2), mar=c(2,2,2,2))
for (i in 1:ncol(X)) {
    param_values = lapply(samples, function(l) l[[i]])
    param_values = unlist(param_values)

    hist(
        param_values,
        ylim=c(0,30000),
        xlim=range(mleest[i], param_values),
        main=names(mhout6$beta_mat[[1]])[i]
    )
    abline(v=mleest[i], col="blue", lwd=2, lty=2)
}

