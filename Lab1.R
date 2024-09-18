###############
# Q1 a)

n = 65
sumx = 24890
alpha = 1
beta = 0.01

pmean = (alpha+sumx)/(beta+n)

L = qgamma(0.025, shape = alpha+sumx, rate = beta+n)
U = qgamma(0.975, shape = alpha+sumx, rate = beta+n)

cat("Posterior mean: ",pmean," (",L,",",U,")")

###############
# Q1 b)

N = 500
L = rep(NA, length=250)
U = rep(NA, length=250)

for (i in 1:250) {
  data = sort(rgamma(N, shape = alpha+sumx, rate = beta+n))
  L[i] = data[0.025*N]
  U[i] = data[0.975*N]
}

widthL=max(L)-min(L)
widthU=max(U)-min(U)
par(mfrow=c(1,2))
hist(L,probability=T,xlab="Lower 95% CI bound")
points(qgamma(0.025, shape = alpha+sumx, rate = beta+n),0,pch=16,col=2)
hist(U,probability=T,xlab="Upper 95% CI bound")
points(qgamma(0.975, shape = alpha+sumx, rate = beta+n),0,pch=16,col=2)
cat("L interval variability (range):",widthL,"\n")
cat("U interval variability (range):",widthU,"\n")

###############
# Q1 c)

n = 65
sumx = 24890
alpha = 1
beta = 0.01

theta = rgamma(1000, shape = alpha+sumx, rate = beta+n)
y = rpois(1000, theta)

hist(y, probability=T, ylab="Density", main="Posterior Predictive Distribution")

xvals = 300:450

predDist=dnbinom(xvals, alpha+sumx, 1-1/(beta + n + 1))

for (i in 1:length(xvals)) {
  lines(xvals, predDist, col="purple", lwd=2)
}

###############
# Q1 d) in onenote



###############
# Q3 a)
# set d = 1 and ell = 1 for example

buffon=function(n, d=1, ell=1, make.plot=TRUE) {
  #### runif is uniform random sampling
  x=runif(n,0,d) # generate samples of x
  theta=runif(n,0,pi) # generate samples of theta
  
  #### issue with experiment, need pi to generate samples of theta but we're
  #### trying to estimate pi so it's a bit backward
  
  
  L=sin(theta)
  int=(ell*L)>x # determine samplesd x values with correspond to an intersection
  C=cumsum(int) # calculate counts of intersection for each value of n
  
  if (make.plot) {
    plot(1:n, (1:n)*2*ell/(d*C), type="l", xlab="Number of simulations",
         ylab=expression(paste("Estimate of ", pi)), ylim=c(2,6))
    abline(pi, 0, col=2)
  }
  return(n*2*ell/(d*C[n]))
}

buffon(2000)

###############
# Q3 b) in onenote

# Q3 c)

out = rep(NA, 1000)
n=2000
ell.seq=seq(0.1, 1, length=20)
sd.out=rep(NA, length(ell.seq))
for (j in 1:length(ell.seq)) {
  for (i in 1:1000) {
    out[i]=buffon(n, ell=ell.seq[j], make.plot = F)
    sd.out[j]=sd(out)
  }
}
plot(ell.seq, sd.out, xlab = "l", ylab = "Estimate Standard Error", type="l")
