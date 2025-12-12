rm(list = ls())  # Clear environment
gc()             # Free memory

# setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")

library(MASS)
library(latex2exp)
library(squash)
library(DRPT)

set.seed(110932)

# Set parameters
n = 250
m = 250
d = 1

MC = 300
S = 50

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

invCDF = function(x){
  return(sqrt(x))
}


#---------------------------------
# Estimate LL
#---------------------------------
size_n = numeric(4)
n_index = 1
alpha = 0

for (n_train in c(100, 200, 1000, 2000)){
  # r_N  = r = function(x,y) {
  #   return(2*x*(1+50*sin(50*x)/n_train))
  # }
  
  r_N  = r = function(x,y) {
    return(2*x^(1+500/n_train))
  }
  x = runif(1000)
  plot(x, r_N(x,x))
}

for (n_train in c(100, 200, 1000, 2000)){
  # r_N  = r = function(x,y) {
  #   return(2*x*(1+50*sin(50*x)/n_train))
  # }
  
  r_N  = r = function(x,y) {
    return(2*x^(1+100/n_train))
  }
  
  sum_ind = 0
  
  #--------------------------------------
  # Monte Carlo simulation on test data
  #--------------------------------------
  for (b in 1:MC) {
    ######## generate data
    X = runif(n, 0, 1)
    Y = invCDF(runif(m, 0, 1))
    
    ########### Our test: call the SPT function
    p_M_alpha = DRPT(X, Y, r = r_N, kernel = gaussian.kernel)
    
    if (p_M_alpha < 0.05) {
      sum_ind = sum_ind + 1
    }
  }
  
  # Calculate the power for this alpha
  size_n[n_index] = sum_ind / MC
  print(sum_ind / MC)
  n_index = n_index + 1
}

# Save the data to a CSV file
data_to_save = data.frame(n = c(100, 200, 1000, 2000), Size = size_n)
write.csv(data_to_save, "experiments/results/r_estimation_r_N.csv", row.names = FALSE)