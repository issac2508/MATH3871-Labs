#Template Script for Assignment, MATH3871/5960 

#install.packages("mvtnorm")
library(MASS) #you may need to install the package first! 
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

wine = read.csv("winequality-red.csv")

anyNA(wine)
# returned false so there are no NA values

#-------------------------------------------------------------------------------
#QUESTION 2.

wine$good = as.integer(wine$quality >= 6.5)
wine = subset(wine, select = -quality)

#-------------------------------------------------------------------------------
#QUESTION 3. 

model = glm(good ~ ., data=wine, family=binomial)
mleest = as.numeric(coef(model))

#-------------------------------------------------------------------------------

############################# PART II Bayesian #################################
#-------------------------------------------------------------------------------
#QUESTION 2.

lpost.LR <- function(beta,X,y) {
  k = ncol(X)
  alpha = numeric(k)
  omega = 100 * diag(k)
  
  
  loglikelihood = y%*%X%*%beta-as.matrix(sum(log(1+exp(X%*%beta))))
  # find logprior
  logprior = -0.5 * t(beta - alpha) %*% solve(omega) %*% (beta - alpha) -
    (k / 2) * log(2 * pi) - 0.5 * log(det(omega))
  
  logpost = loglikelihood + logprior
  
  return(logpost)
}

#-------------------------------------------------------------------------------
#QUESTION 3. 

mhmcmc <- function(y, X, B, nsims, Sigma) {
  
  k = ncol(X)
  beta_mat = matrix(0, nrow = nsims, ncol = k)
  accprob = numeric(nsims)
  beta_star = matrix(0, nrow = nsims, ncol = k)
  acc = 0
  
  beta = B
  
  for (i in 1:nsims) {
    prop_beta = mvrnorm(1, mu=beta, Sigma=Sigma)
    
    log_post_current = lpost.LR(beta, X, y)
    log_post_proposed = lpost.LR(prop_beta, X, y)

    accprob[i] = exp(log_post_proposed - log_post_current)
    if (is.na(accprob[i])) {
      acceptance_threshold = 1
    } else {
      acceptance_threshold = min(1, accprob[i])
      
    }
    
    if (runif(1) <= acceptance_threshold) {
      beta = prop_beta  # Accept the proposed value
      acc = acc + 1     # Increment accepted count
    }
    
    beta_mat[i, ] = beta
    beta_star[i, ] = prop_beta
  }
  
  #return variables 
  mhout = list("beta_mat"=beta_mat, "accprob"=accprob, "beta_star" = beta_star, "acc"=acc)
  return(mhout)
}

#-------------------------------------------------------------------------------
#QUESTION 4. 

#Covariance for proposal: 
Sigma = diag(c(1100, 0.0015, 0.04, 0.05, 0.0005, 0.4, 0.00001, 0.000001, 0.0001, 0.09, 0.03, 0.002)) #Do not modify  

X=cbind(1, as.matrix(subset(wine, select = -good)))
y=t(as.matrix(wine$good))

mhout1 = mhmcmc(y=y, X=X, B=mleest, nsims=10^5, Sigma=Sigma)

#(a) produce trace plots:

num_coefficients = ncol(mhout1$beta_mat)
names = c(list("intercept"), colnames(subset(wine, select = -good)))

par(mfrow=c(4,3))
for (i in 1:num_coefficients) {
  plot(mhout1$beta_mat[, i], type = "l", main = paste("Trace Plot for", names[i]),
       xlab = "Iteration", ylab = "Coefficient Value")
}

#(b) proportion of accepted moves:
mh4acc = mhout1$acc/10^5

#-------------------------------------------------------------------------------
#QUESTION 5.

init1 = mleest + rnorm(length(mleest), 0, 10)
init2 = matrix(50, nrow=12, ncol=1)
init3 = matrix(-50, nrow=12, ncol=1)
init4 = matrix(0, nrow=12, ncol=1)

chain1 = mhmcmc(y=y, X=X, B=init1, nsims=10^5, Sigma=Sigma)
chain2 = mhmcmc(y=y, X=X, B=init2, nsims=10^5, Sigma=Sigma)
chain3 = mhmcmc(y=y, X=X, B=init3, nsims=10^5, Sigma=Sigma)
chain4 = mhmcmc(y=y, X=X, B=init4, nsims=10^5, Sigma=Sigma)

par(mfrow = c(4, 3), oma = c(0, 0, 0, 6))
for (i in 1:num_coefficients) {
  plot(chain1$beta_mat[, i], type = "l", col = "blue", lty = 1,
       main = paste("Trace Plot for", names[i]),
       xlab = "Iteration", ylab = "Coefficient Value",
       ylim = range(c(chain1$beta_mat[, i], chain2$beta_mat[, i],
                      chain3$beta_mat[, i], chain4$beta_mat[, i])))
  
  lines(chain2$beta_mat[, i], col = "red", lty = 1)
  lines(chain3$beta_mat[, i], col = "green", lty = 1)
  lines(chain4$beta_mat[, i], col = "purple", lty = 1)
  
}
par(xpd = NA)
legend("topright", inset = c(-0.85, -20), legend = c("Chain 1", "Chain 2", "Chain 3", "Chain 4"),
       col = c("blue", "red", "green", "purple"), lty = c(1, 1, 1, 1), cex = 0.8)


#-------------------------------------------------------------------------------
#QUESTION 6.

cov_matrix = as.matrix(vcov(model))
cov_mat = 0.05 * cov_matrix

mhout2 = mhmcmc(y=y, X=X, B=mleest, nsims=10^5, Sigma=cov_mat)

#calculate the acceptance rate
mh6acc = mhout2$acc/10^5

#(a)produce trace plots:
par(mfrow=c(4,3), oma = c(0, 0, 0, 0))
for (i in 1:num_coefficients) {
  plot(mhout2$beta_mat[, i], type = "l", main = paste("Trace Plot for", names[i]),
       xlab = "Iteration", ylab = "Coefficient Value")
}

#(b) produce marginal histograms and overlay MLE: 
par(mfrow=c(4,3))
for (i in 1:num_coefficients) {
  acf(mhout2$beta_mat[, i], main = paste("ACF for", names[i]))
}

post_burnin = mhout2$beta_mat[-(1:100), ]

par(mfrow=c(4,3))
for (i in 1:num_coefficients) {
  x_limits = range(post_burnin[, i], mleest[i])
  
  hist_data = hist(post_burnin[, i], breaks = 50, main = paste("Histogram for", names[i]),
       xlab = "Coefficient Value", col = "lightblue", probability = TRUE, xlim = x_limits)
  
  segments(x0 = mleest[i], y0 = 0, x1 = mleest[i], y1 = max(hist_data$density), 
           col = "red", lwd = 1, lty = 1)
}

