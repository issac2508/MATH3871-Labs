###############
# Q3 a)
# set d = 1 and ell = 1 for example

buffon=function(n, d=1, ell=1, make.plot=TRUE) {
  #### runif is uniform random sampling
  x=runif(n,0,d) # generate samples of x
  theta=runif(n,0,pi) # generate samples of theta
  
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
# Q3 b)

