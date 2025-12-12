rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")

library(latex2exp)
library(DRPT)


#------------------------------------------------------------------------------
# d = 4 - U/V/D
#------------------------------------------------------------------------------
n = 250
m = 250
r = c(1, 1, 1, 1)
r.func = function(x){
  return(r[x+1])
}
MC = 5000

# Preallocate power and standard deviation vectors for efficiency
pow_vector_U = numeric(15)
pow_vector_V = numeric(15)
pow_vector_D = numeric(15)

# Precompute the alpha values
alpha_values = seq(from = 0, to = 0.4, length.out = 15)

# Main loop over alpha values
for (i in seq_along(alpha_values)) {
  alpha = alpha_values[i]
  gamma = (25*alpha/43) * (1/sqrt(3) + 1/sqrt(10))
  print(gamma)
  
  # Preallocate the power for each iteration of MC
  pow_alpha_U = numeric(MC)
  pow_alpha_V = numeric(MC)
  pow_alpha_D = numeric(MC)
  
  
  for (M in 1:MC) {
    # Sampling for X and Y based on the probabilities
    X = sample(0:3, n, prob = c(1/4, 1/4, 1/4, 1/4), replace = TRUE)
    Y = sample(0:3, m, prob = c(1/4, 1/4, (1+gamma)/4, (1-gamma)/4), 
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
data_to_save =data.frame(Alpha = alpha_values, Power_V = pow_vector_V,  
                           Power_U = pow_vector_U, Power_D = pow_vector_D)
write.csv(data_to_save, "experiments/results/4d_1111_discrete_data_UVD.csv", row.names = FALSE)

data.discrete_UV = read.csv(file = "experiments/results/4d_1111_discrete_data_UVD.csv")

png("experiments/pictures/4d_1111Appendix.png")
par(mfrow = c(1, 1))
plot(data.discrete_UV$Alpha, data.discrete_UV$Power_U, type = "b", col = "darkviolet", lwd = 2,
     pch = 14, xlab = TeX(r'($\eta$)', bold = TRUE),
     ylab = "Power", ylim = c(0,1))
lines(data.discrete_UV$Alpha, data.discrete_UV$Power_V, type = "b", col = "blue", lwd = 2,
      pch = 15)
lines(data.discrete_UV$Alpha, data.discrete_UV$Power_D, type = "b", col = "cyan", lwd = 2,
      pch = 15)
legend("topleft", inset = +0.01,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("Discrete DRPT - U", "Discrete DRPT - V", "Discrete DRPT - D"),
       col = c("darkviolet", "blue", "cyan"),
       pch = c(14, 14, 15))
abline(h = 0.05, col = "red", lty = 2)
dev.off()


 #### BOTH PLOTS TOGETHER #########
# --- Read data ---
# Per your note:
#   1) 4d_1111_discrete_data_UVD.csv -> r = (1,3,3,10)
#   2) 4d_13310_discrete_data_UVD.csv -> r = (1,1,1,1)
data_r13310 =read.csv("experiments/results/4d_1111_discrete_data_UVD.csv")  # r = (1,3,3,10)
data_r1111  =read.csv("experiments/results/4d_13310_discrete_data_UVD.csv") # r = (1,1,1,1)

# --- Output file ---
png("experiments/pictures/4d_combinedAppendix.png", width = 1200, height = 900, res = 150)

par(mfrow = c(1,1), mar = c(5,5,2,2))

# Colors by statistic
colV ="darkviolet"; colU ="blue"; colD ="cyan"

# Plot base using r = (1,1,1,1), linetype = solid
plot(data_r1111$Alpha, data_r1111$Power_V, type = "b",
     col = colU, lwd = 2, pch = 16, lty = 1,
     xlab = TeX(r'($\eta$)', bold = TRUE), ylab = "Power",
     ylim = c(0, 1))
lines(data_r1111$Alpha, data_r1111$Power_U, type = "b",
      col = colV, lwd = 2, pch = 16, lty = 1)
lines(data_r1111$Alpha, data_r1111$Power_D, type = "b",
      col = colD, lwd = 2, pch = 16, lty = 1)

# Add r = (1,3,3,10), linetype = dashed
lines(data_r13310$Alpha, data_r13310$Power_V, type = "b",
      col = colU, lwd = 2, pch = 1, lty = 2)
lines(data_r13310$Alpha, data_r13310$Power_U, type = "b",
      col = colV, lwd = 2, pch = 1, lty = 2)
lines(data_r13310$Alpha, data_r13310$Power_D, type = "b",
      col = colD, lwd = 2, pch = 0, lty = 2)

# Reference line
abline(h = 0.05, col = "red", lty = 2)

# Legends: one for statistics (colors), one for r (linetypes)
legend("topleft", inset = 0.01, bty = "n", lty = 1, lwd = 2, horiz = FALSE,
       legend = c("(E4)", "(E5)", "(E6)"),
       col = c(colV, colU, colD), pch = c(16,16,16), title = "Statistic (color)")

legend("left", inset = 0.02, bty = "n", lty = c(1,2), lwd = 2, horiz = FALSE,
       legend = c("r = (1,1,1,1)", "r = (1,3,3,10)"),
       title = "r (linetype)")

dev.off()
