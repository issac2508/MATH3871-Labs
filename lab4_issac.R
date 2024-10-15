##########
# Q1

# a)
gx=runif(1000)
w=gx^2*exp(-gx) # weights
nw=w/sum(w) # normalised weights

# b)
lines(density(gx,weights=nw),col=4,lty=2)


##########
# Q3

# a)
xx=seq(0,1,length=150)
plot(xx, dnorm(xx,0.5,0.01), type="l", xlab="u", ylab="density")
u=runif(1000)
w=dnorm(u, 0.5, 0.01)
ess=function(w) {w=w/sum(w); return(1/sum(w^2))}
cat("The ESS is", ess(w), "\n")

# b)
L = 100
nvec = seq(1000, 200000, length=L)
ESS=rep(NA, L)

for (i in 1:length(nvec)) {
  u=runif(nvec[i])
  ESS[i]=ess(dnorm(u, 0.5, 0.01))
}

index = which(ESS>5000)[1]
cat("Number of samples is", nvec[index])

# c)
## Step 1
p=0.02
u=runif(1000)
w=dnorm(u, 0.5, 0.01)^p
cat("The ESS is", ess(w), "\n")

xx=seq(0,1,length=250)
a=sample(u, 10000, replace=T, prob=w/sum(w)) # generated equally weighted sample
hist(a, xlab="u", ylab="density", probability=T, ylim=c(0, dnorm(0,0,0.01)/5), 
     main="Sequential Importance Sampling")
lines(xx, dnorm(xx, 0.5, 0.01)/5, type="l")
lines(xx, dnorm(xx, 0.5, 0.01)^p*5, lty=2, col=2)
m.a=mean(a)
sd.a=sd(a)

## Step 2
x=rnorm(1000, m.a, sd.a) # new proposal distribution g2(x)
w2=dnorm(x,0.5,0.01)/dnorm(x,m.a,sd.a)
cat("The ESS is", ess(w2), "\n")

# d)
L = 100
nvec = seq(1000, 200000, length=L)
ESS=rep(NA, L)

for (i in 1:length(nvec)) {
  x=rnorm(nvec[i], m.a, sd.a) # new proposal distribution g2(x)
  w2=dnorm(x,0.5,0.01)/dnorm(x,m.a,sd.a)
  ESS[i]=ess(w2)
}

index = which(ESS>5000)[1]
cat("Number of samples is", nvec[index] + 1000)

# e)
# since the first method needed ~140,000 samples and the second
# method only needs ~30,000 samples, the second procedure is much
# more effective

##########
# Q4