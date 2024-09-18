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
# Q3 b) in one note

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
