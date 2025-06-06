rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(MASS)
library(latex2exp)
library(squash)

source(file = "DRPT.R")

##################################################################
# set parameters
##################################################################

n = 250
m = 250
d = 2

MC = 500
S = 50
xxx = 15

r = function(x,y) {
  return(2*x)
}

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

invCDF = function(x){
  return(sqrt(x))
}

##################################################################
# Simulation for varying alpha
##################################################################
pow_vector = numeric(xxx)
alpha_index = 1

for (alpha in seq(from = 0, to = 0.9, length.out = xxx)){
  print(alpha)
  sum_ind = 0
  
  for (b in 1:MC) {
    X = cbind(runif(n, 0, 1), runif(n, 0, 1))
    
    Y = matrix(0, nrow = m, ncol = 2)
    binary_vec = rbinom(m, 1, 1 / (1 + alpha))
    Y[binary_vec == 1, ] = cbind(invCDF(runif(sum(binary_vec))), runif(sum(binary_vec)))
    Y[binary_vec == 0, ] = cbind(rbeta(sum(!binary_vec), 0.5, 0.5), rbeta(sum(!binary_vec), 0.5, 0.5))
    
    median.bandwidth = median(as.vector(dist(rbind(X,Y),
                                             method = "euclidean", diag = FALSE, upper = FALSE)))
    
    # Our test: call the DRPT function
    p_M_alpha = DRPT(X, Y, S, H = 99, r, 
                    kernel = function(x,y) gaussian.kernel(x, y, lambda = median.bandwidth))
    
    if (p_M_alpha < 0.05) {
      sum_ind = sum_ind + 1
    }
  }
  
  pow_vector[alpha_index] = sum_ind / MC
  alpha_index = alpha_index + 1
}

#####################################################################
# Save the data to a CSV file
#####################################################################
data_to_save = data.frame(Alpha = seq(from = 0, to = 0.9, length.out = xxx), 
                           Power = pow_vector)
write.csv(data_to_save, "experiments/results/BIV2simul_shiftedMMD.csv", row.names = FALSE)
