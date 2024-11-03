#################
# MATH3871 Lab 5
#################

# Question 1

n = 50
mu = 5
sig = 2
x_data = rnorm(n, 5, sig)

# X and Y will store the Gibbs samples
N = 1000
X = Y = rep(NA, N)
X[1] = 0.5    # sig^-2 samples
Y[1] = 4.8    # mu samples

alpha = 0.5
beta = 1
xi = 4.7
kappa = 1

# Gibbs sampling
for (i in c(2:N)) {
    X[i] = rgamma(1, alpha+n/2, rate = beta+0.5*sum((x_data-Y[i-1])^2))
    Y[i]=rnorm(1, (X[i]*sum(x_data) + kappa*xi)/(X[i]*n + kappa), sqrt(1/(X[i]*n + kappa)))
}

# Histogram and trace plots
par(mfrow=c(2,2))
hist(X, probability=T)
plot(X, type="l")
hist(Y, probability=T)
plot(Y, type="l")

# Plot samples with 2D contour plot
f=function(x,y) {
    x^(alpha+n/2-1)*exp(-beta*x - (kappa*(y-xi)^2)/2 - (x*sum((x_data - y)^2))/2)
}
xx=seq(0,1,length=100)
yy=seq(3,6,length=100)
zz=matrix(NA,ncol=100,nrow=100)
for (i in 1:100) zz[,i]=f(xx,yy[i])
par(mfrow=c(1,1))
contour(xx,yy,zz,xlab="x",ylab="y",col=2)
points(X,Y,pch=16,cex=0.3)


# Question 2

# a)
# First, we generate 50 samples for x
n = 50
rho = 0.95
b0 = 1
b1 = 1
phi = 10

library(MASS)

mu = rep(0, n)                          # mean
sigma = matrix(rep(rho, n*n), ncol = n) # covariance
diag(sigma) = 1
X = mvrnorm(1, mu, sigma) # num samples is 1 because this 1 sample is an n dimensional vector already
y = mvrnorm(1, b0 + b1*X, 1/phi * diag(n))

# Now, we use the samples to generate the posterior distribution of
# b0, b1, and phi

n0 = 4
S0 = 0.5

xbar = mean(X)
Sxx = sum(X^2) - n * (xbar^2)
Xmat = cbind(rep(1,n), X)   # design matrix: X values with intercept column
bhat = solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% y # `solve` finds the inverse

# Do Gibbs sampling
b0_samples = c()
b1_samples = c()
phi_samples = c()

b0_samples[1] = 0.5 # remember that MCMC initialisation can be any value
b1_samples[1] = 0.5
phi_samples[1] = 1

for (i in c(2:1000)) {
    b0_samples[i] = rnorm(
        1,
        mean = bhat[1] - xbar*(b1_samples[i-1]-bhat[2]),
        sd = sqrt(1/(phi_samples[i-1]*n*Sxx) * (1-n*(xbar^2)/(sum(X^2))))
    )

    b1_samples[i] = rnorm(
        1,
        mean = bhat[2] - xbar*n/sum(X^2) * (b0_samples[i] - bhat[1]),
        sd = sqrt(1/(phi_samples[i-1]*Sxx) * (1 - n*(xbar^2)/sum(X^2)))
    )

    # Sb (S_beta) was computed in tutorial 5 question 5b (apparently)
    Sb = t(y-(Xmat%*%c(b0_samples[i], b1_samples[i])))%*%(y-(Xmat%*%c(b0_samples[i], b1_samples[i])))+n0*S0

    phi_samples[i] = rgamma(1, (n+n0)/2, Sb/2)
}

par(mfrow = c(2,2))
plot(b0_samples, type="l")
plot(b1_samples, type="l")
plot(phi_samples, type="l")

# b)

# Plot autocorrelation of the samples
# - autocorrelation is dependence of each sample on its previous samples
#   (high = samples are highly dependent)

par(mfrow = c(2,2))
acf(b0_samples) # acf = auto correlation function
acf(b1_samples)
acf(phi_samples)

# Determine thinning threshold by where the ACF curves reach 0

# c)

# Reparameterisation is equivalent to centering covariates X

Xnew = X - mean(X)
xbar = mean(Xnew)
Sxx = sum(Xnew^2) - n*(xbar^2)
Xmat = cbind(rep(1,n), Xnew)

ahat = solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% y

a0_samples = c()
a1_samples = c()
phi_samples = c()

a0_samples[1] = 0.5
a1_samples[1] = 0.5
phi_samples[1] = 1

for (i in c(2:1000)) {
    a0_samples[i] = rnorm(
        1,
        mean = ahat[1] - xbar*(a1_samples[i-1]-ahat[2]),
        sd = sqrt(1/(phi_samples[i-1]*n*Sxx) * (1-n*(xbar^2)/(sum(X^2))))
    )

    a1_samples[i] = rnorm(
        1,
        mean = ahat[2] - xbar*n/sum(X^2) * (a0_samples[i] - ahat[1]),
        sd = sqrt(1/(phi_samples[i-1]*Sxx) * (1 - n*(xbar^2)/sum(X^2)))
    )

    Sa = t(y-(Xmat%*%c(a0_samples[i], a1_samples[i])))%*%(y-(Xmat%*%c(a0_samples[i], a1_samples[i])))+n0*S0

    phi_samples[i] = rgamma(1, (n+n0)/2, Sa/2)
}

par(mfrow=c(2,2))
plot(a0_samples, type="l")
plot(a1_samples, type="l")
plot(phi_samples, type="l")

par(mfrow=c(2,2))
acf(a0_samples)
acf(a1_samples)
acf(phi_samples)

# there is basically no autocorrelation here! the reparameterisation provided
# independent components. this means convergence happens faster - more efficient
