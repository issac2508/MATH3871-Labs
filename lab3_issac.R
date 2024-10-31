###########
# 1)

samp=function(n) {
  out=NULL
  K=exp(-1)
  while(length(out)<n) {
    x=runif(1)
    y=runif(1, 0, K)
    ind=(y<x^2*exp(-x))
    out=c(out, x[ind])
  }
  return(out)
}

out=samp(5000)
hist(out, probability=T, xlab="x", ylab="Density")
xx=seq(0,1,length=100)
lines(xx, dgamma(xx,3,1)/pgamma(1,3,1), col="purple")
lines(c(0,1),c(1,1),col="red")

###########
# 2a)
u = runif(10000)

mean(exp(exp(u)))

###########
# 2b)

2*mean(exp(-(1/u-1)^2)/(u^2))

###########
# 2b)

###########
# 2d)
u1 = runif(10000)
u2 = runif(10000)

mean(1/u1^2*1/u2^2*exp(2-1/u1-1/u2)*(u1<u2))

###########
# 3)

L=20
nvec=seq(100, 20000, length=L)
mvec=vvec=rep(NA, L)
for (i in 1:L) {
  u=runif(nvec[i])
  hu=(cos(50*u)+sin(20*u))^2
  mvec[i]=mean(hu)
  vvec[i]=1/(nvec[i])^2 * sum((hu-mvec[i])^2)  
}

plot(nvec, mvec, type="l", ylim=c(0.9, 1.1), xlab="m", ylab="Estimate")
lines(nvec, mvec+sqrt(vvec), lty=2, col="purple")
lines(nvec, mvec-sqrt(vvec), lty=2, col="purple")
abline(h=0.965, col="red")



