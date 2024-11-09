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
  
  beta = matrix(0, nrow=k, ncol=1)
  
  
  for (i in 1:nsims) {
    prop_beta = mvrnorm(1, mu=beta, Sigma=Sigma)

    log_post_current = lpost.LR(beta, X, y)
    log_post_proposed = lpost.LR(prop_beta, X, y)
    
    print(log_post_proposed)
    print(log_post_current)

    accprob[i] = exp(log_post_proposed - log_post_current)
    acceptance_threshold = min(1, accprob[i])
    
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



#-------------------------------------------------------------------------------
#QUESTION 6. 



#calculate the acceptance rate
mh6acc = 

#(a)produce trace plots:


#(b) produce marginal histgrams and overlay MLE: 


