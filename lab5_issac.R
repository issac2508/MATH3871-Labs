
##########
# Q2

# a)
# params
rho=0.95
phi=10
n0=4
s0=0.5
beta0=1
beta1=1

# simulate x

library(MASS)
n=50
Sigma=matrix(rep(rho, n*n), ncol=n)
diag(Sigma)=1

X=mvrnorm(1, rep(0,n), Sigma)
Y=mvrnorm(1, beta0 + beta1*X, 1/phi*diag(n))

xbar=mean(X)
Sxx=sum(X^2)-n*(xbar^2)
Xmat=cbind(rep(1,n), X)
betahat=solve(t(Xmat)%*%Xmat)%*%t(Xmat)%*%Y

niter=1000
beta0MCMC=c()
beta1MCMC=c()
phiMCMC=c()
beta0MCMC[1]=0.5 # some initial values
beta1MCMC[1]=0.5
phiMCMC[1]=1


# gibbs sampling
for (i in c(2:niter)) {
  beta0MCMC[i]=rnorm(1, betahat[1]-xbar*(beta1MCMC[i-1]-betahat[2]), sqrt(1/phiMCMC[i-1]/Sxx/n*(1-n*(xbar^2)/sum(X^2))))
  beta1MCMC[i]=rnorm(1, betahat[2]-xbar*n/sum(X^2)*(beta0MCMC[i]-betahat[1]), sqrt(1/phiMCMC[i-1]/Sxx*(1-n*(xbar^2)/sum(X^2)))) 
  
  Sbeta=t(Y-(Xmat%*%c(beta0MCMC[i], beta1MCMC[i])))%*%(Y-(Xmat%*%c(beta0MCMC[i], beta1MCMC[i])))+n0*s0
  phiMCMC[i]=rgamma(1, ((n+n0)/2), (Sbeta/2))
  
}

par(mfrow=c(2,2))
plot(beta0MCMC, type="l")
plot(beta1MCMC, type="l")
plot(phiMCMC, type="l")

# b)
par(mfrow=c(2,2))
acf(beta0MCMC)
acf(beta1MCMC)
acf(phiMCMC)

# c)
Xnew=X-mean(X)
xbar=mean(Xnew)
Sxx=sum(Xnew^2)-n*(xbar^2)
Xmat=cbind(rep(1,n), Xnew)
alphahat=solve(t(Xmat)%*%Xmat)%*%t(Xmat)%*%Y

niter=1000
alpha0MCMC=c()
alpha1MCMC=c()
phiMCMC=c()
alpha0MCMC[1]=0.5
alpha1MCMC[1]=0.5
phiMCMC[1]=1

for (i in c(2:niter)) {
  alpha0MCMC[i]=rnorm(1, alphahat[1]-xbar*(alpha1MCMC[i-1]-alphahat[2]), sqrt(1/phiMCMC[i-1]/Sxx/n*(1-n*(xbar^2)/sum(Xnew^2))))
  alpha1MCMC[i]=rnorm(1, alphahat[2]-xbar*n/sum(Xnew^2)*(alpha0MCMC[i]-alphahat[1]),sqrt(1/phiMCMC[i-1]/Sxx*(1-n*(xbar^2)/sum(Xnew^2))))
  
  Sbeta=t(Y-(Xmat%*%c(alpha0MCMC[i], alpha1MCMC[i])))%*%(Y-(Xmat%*%c(alpha0MCMC[i], alpha1MCMC[i])))+n0*s0
  
  phiMCMC[i]=rgamma(1,((n+n0)/2), (Sbeta/2))
}

beta1MCMC=alpha1MCMC
beta0MCMC=alpha0MCMC-mean(X)*beta1MCMC
par(mfrow=c(2,2))
plot(beta0MCMC, type="l")
plot(beta1MCMC, type="l")
plot(phiMCMC, type="l")

par(mfrow=c(2,2))
acf(beta0MCMC)
acf(beta1MCMC)
acf(phiMCMC)
