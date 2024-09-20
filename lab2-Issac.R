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

obs = inverse_func(runif(5000)) # do in array form

hist(obs, probability=T)
