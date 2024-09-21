##########
# Q1)

inverse_func = function(x) {
  return ((-1+sqrt(1+8*x))/2)
}

obs = inverse_func(runif(5000))

hist(obs, probability=T)

xvals=seq(0,1,length=100)

lines(xvals, xvals+0.5, col="purple")

###########
# Q2)

inverse_func = function(x) {
  if (x < 0.5) {
    return ((log(2*x, base=exp(1)))/2)
  } else {
    return (-(log(2*(1-x), base=exp(1)))/2)
  }
}

obs = rep(NA, length=5000)
samples = runif(5000)

for (i in 1:5000) {
  obs[i] = inverse_func(samples[i])
}

hist(obs, probability=T, breaks = 40, xlim=c(-4,4), ylim=c(0,1))
# breaks changes the number of bins
# xlim and ylim change the scale of the graph

xvals=seq(-5,5,length=100)
lines(xvals, exp(-sign(xvals)*2*xvals), col="purple")

###########
# Q3)
gen_samples = function(num, lambda) {
  return(-log(1-runif(num))/lambda)
}

lambda = 1 # choose random?

par(mfrow=c(3,3), oma = c(0, 0, 2, 0))

nums=c(10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000)
xvals=seq(0,10,length=100)

for (num in nums) {
  hist(gen_samples(num, lambda), probability=T, breaks = 20, xlim=c(0,10),
       ylim=c(0,1), main = paste("num samples = ", as.character(num)))
  lines(xvals, lambda * exp(-lambda*xvals), col="purple")
}

mtext(paste("Lambda = ", as.character(lambda)), outer = TRUE)

###########
# Q4)



