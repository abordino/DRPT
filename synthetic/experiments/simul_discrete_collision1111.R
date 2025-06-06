rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(latex2exp)

source(file = "DRPT.R")

#####################################################################
# setting the parameters for d = 4, (1,1,1,1), MMD
#####################################################################
n = 250
m = 250

v = c(1, 3, 3, 10)
r = function(i){
  return(v[i+1])
}

collision.kernel = function(x, y){
  return(sum((x == 0:3)*(y == 0:3)))
}

##################################################################
# Simulation for varying alpha
##################################################################

MC = 500; S = 50
pow_vector = numeric(30)
alpha_values = seq(from = 0, to = 0.8, length.out = 30)

for (i in seq_along(alpha_values)) {
  alpha = alpha_values[i]
  gamma = (25*alpha/43) * (1/sqrt(3) + 1/sqrt(10))
  print(alpha)
  
  pow_alpha = numeric(MC)
  
  for (b in 1:MC) {
    
    # Sampling for X and Y based on the probabilities
    X = sample(0:3, n, prob = c(1/4, 1/4, 1/4, 1/4), replace = TRUE)
    Y = sample(0:3, m, prob = c(1/4, 1/4, (1+gamma)/4, (1-gamma)/4), 
               replace = TRUE)
    
    # Compute p-value using DRPT
    p_M_alpha = DRPT(X, Y, S, H = 99, r, kernel = collision.kernel)
    pow_alpha[b] = (p_M_alpha < 0.05)
    
  }

  pow_vector[i] = mean(pow_alpha)
}

#####################################################################
# Save the data to a CSV file
#####################################################################
data_to_save = data.frame(Alpha = alpha_values, Power = pow_vector)
write.csv(data_to_save, "experiments/results/4d_1111_MMD.csv", row.names = FALSE)

