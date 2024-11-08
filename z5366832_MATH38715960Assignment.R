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
    prior = dmvnorm(x=beta, mean=alpha, sigma=omega)

    log.likelihood = t(y) %*% X %*% beta - sum(log(1 + exp(X %*% beta)))

    return (log(prior) + log.likelihood)
}

#### Just making sure it works
y = as.matrix(wine$good)
X = wine[,!names(wine) %in% c("quality", "good")]
X = as.matrix(cbind("intercept"=rep(1, nrow(X)), X))
beta = mleest

lpost.LR(beta, X, y)


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
        # if (lpost.LR(B, X, y) > 0) {
        #     alpha = min(1, lpost.LR(B_prop, X, y)/lpost.LR(B, X, y))
        # } else {
        #     alpha = 1
        # }
        # log.alpha = log(alpha)

        log.alpha = lpost.LR(B_prop, X, y) - lpost.LR(B, X, y)
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

nSims = 10^5 # change to 10^5

mhout1 = mhmcmc(y, X, B, nSims, Sigma)

#(a) produce trace plots:
mhout1

# Just doing for a single i to see how to access all the values
# i = 1
# param_values = lapply(mhout1$beta_mat, function(l) l[[i]])
# param_values = unlist(param_values)
#
# par(mfrow=c(1,1))
# plot(param_values, type="l")
# param_values

par(mfrow=c(6,2))
for (i in 1:ncol(X)) {
    # for each ith coefficient...
    param_values = lapply(mhout1$beta_mat, function(l) l[[i]]) # get list of accepted values for the current coefficient
    param_values = unlist(param_values)

    plot(param_values, type="l", main=names(mhout1$beta_mat[[1]])[i], xlab="iteration")
}

#(b) proportion of accepted moves:
mh4acc = mhout1$acc / nSims
mh4acc

#-------------------------------------------------------------------------------
#QUESTION 5.



#-------------------------------------------------------------------------------
#QUESTION 6.



#calculate the acceptance rate
mh6acc =

#(a)produce trace plots:


#(b) produce marginal histgrams and overlay MLE:


