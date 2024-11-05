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
  
  loglikelihood = y%*%X%*%beta-log(1+exp(x%*%beta))
  # find logprior
  logprior = -0.5 * t(beta - alpha) %*% solve(omega) %*% (beta - alpha) -
    (k / 2) * log(2 * pi) - 0.5 * log(det(omega))
  
  logpost = loglikelihood + logprior
  
  return(logpost)
}			

#-------------------------------------------------------------------------------
#QUESTION 3. 

mhmcmc <- function(y, X, B, nsims, Sigma) {

  #You must use mvrnorm for random sampling from a multivariate normal. 
  
  
}


 
#return variables 
mhout = list("beta_mat"=beta_mat, "accprob"=accprob, "beta_star" = beta_star, "acc"=acc)
return(mhout)

}



#-------------------------------------------------------------------------------
#QUESTION 4. 

#Covariance for proposal: 
Sigma = diag(c(1100, 0.0015, 0.04, 0.05, 0.0005, 0.4, 0.00001, 0.000001, 0.0001, 0.09, 0.03, 0.002)) #Do not modify  

mhout1 = mhmcmc(  )

#(a) produce trace plots:


#(b) proportion of accepted moves:
mh4acc = 


#-------------------------------------------------------------------------------
#QUESTION 5. 



#-------------------------------------------------------------------------------
#QUESTION 6. 



#calculate the acceptance rate
mh6acc = 

#(a)produce trace plots:


#(b) produce marginal histgrams and overlay MLE: 


