rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")

library(latex2exp)
library(DRPT)

#------------------------------------------------------------------------------
# d = 4 - V and U
#------------------------------------------------------------------------------
n = 250
m = 250
r = c(1, 3, 3, 10)
r.func = function(x){
  return(r[x+1])
}
MC = 5000

# Preallocate power and standard deviation vectors for efficiency
pow_vector_V = numeric(15)
pow_vector_U = numeric(15)
pow_vector_D = numeric(15)

# Precompute the alpha values
alpha_values = seq(from = 0, to = 0.4, length.out = 15)

# Main loop over alpha values
for (i in seq_along(alpha_values)) {
  alpha = alpha_values[i]
  print( (25*alpha/43) * (1/sqrt(3) + 1/sqrt(10))/8 )
  
  # Preallocate the power for each iteration of MC
  pow_alpha_V = numeric(MC)
  pow_alpha_U = numeric(MC)
  pow_alpha_D = numeric(MC)
  
  for (M in 1:MC) {
    # Sampling for X and Y based on the probabilities
    X = sample(0:3, n, prob = c(1/8, 1/8, 3/8, 3/8), replace = TRUE)
    Y = sample(0:3, m, prob = c(1/43, 3/43, (9+25*alpha)/43, (30-25*alpha)/43), 
               replace = TRUE)
    
    # Compute p-value using SPT
    p_M_alpha_U = discrete.DRPT(X, Y, r, H = 99, type = "U") 
    p_M_alpha_V = discrete.DRPT(X, Y, r, H = 99, type = "V") 
    p_M_alpha_D = discrete.DRPT(X, Y, r, H = 99, type = "D") 
    
    # Store whether the p-value is below 0.05 (for power calculation)
    pow_alpha_U[M] = (p_M_alpha_U < 0.05)
    pow_alpha_V[M] = (p_M_alpha_V < 0.05)
    pow_alpha_D[M] = (p_M_alpha_D < 0.05)
  }
  
  # Store mean and standard deviation of power for this alpha
  pow_vector_U[i] = mean(pow_alpha_U)
  pow_vector_V[i] = mean(pow_alpha_V)
  pow_vector_D[i] = mean(pow_alpha_D)
}

# Save the data to a CSV file
data_to_save <- data.frame(Alpha = alpha_values, Power_V = pow_vector_V,  
                           Power_U = pow_vector_U, Power_D = pow_vector_D)
write.csv(data_to_save, "experiments/results/4d_13310_discrete_data_UVD.csv", row.names = FALSE)

data.discrete_UV = read.csv(file = "experiments/results/4d_13310_discrete_data_UVD.csv")

png("experiments/pictures/4d_13310Appendix.png")
par(mfrow = c(1, 1))
plot(data.discrete_UV$Alpha, data.discrete_UV$Power_V, type = "b", col = "darkviolet", lwd = 2,
     pch = 14, xlab = TeX(r'($\eta$)', bold = TRUE),
     ylab = "Power")
lines(data.discrete_UV$Alpha, data.discrete_UV$Power_U, type = "b", col = "blue", lwd = 2,
      pch = 14, xlab = TeX(r'($\eta$)', bold = TRUE),
      ylab = "Power")
lines(data.discrete_UV$Alpha, data.discrete_UV$Power_D, type = "b", col = "cyan", lwd = 2,
      pch = 15)
legend("center", inset = +0.04,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("(E3)", "(E4)", 
                  "(E5)"),
       col = c("darkviolet", "blue", "cyan"),
       pch = c(14, 14, 15))
abline(h = 0.05, col = "red", lty = 2)
dev.off()
