rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")

library(latex2exp)
library(rootSolve)
library(DRPT)

set.seed(200801)

#------------------------------------------------------------------------------
# d = 2, various m 
#------------------------------------------------------------------------------
n = 100
MC = 5000

r = c(1, 3)
r.func = function(x){
  return(r[x+1])
}

# Preallocate a list to store data for each `m` value
all_results = list()

# Preparing the plot
png("experiments/pictures/binary.png")
par(mfrow = c(1, 1))
plot(1, type = "n",  xlim = c(0,0.9), ylim = c(0,1),
     xlab = TeX(r'($\eta$)', bold = TRUE), 
     ylab = "Power")

ind = 1
for (m in c(100, 200, 500, 2000)){
# for (m in c(100)){
  # Preallocate power and standard deviation vectors for efficiency
  pow_vector = numeric(30)
  
  # Precompute the alpha values
  alpha_values = seq(from = 0, to = 0.9, length.out = 30)
  
  # Main loop over alpha values
  for (i in seq_along(alpha_values)) {
    alpha = alpha_values[i]
    
    # Preallocate the power for each iteration of MC
    pow_alpha = numeric(MC)
    
    for (M in 1:MC) {
      # Sampling for X and Y based on the probabilities
      X = rbinom(n, 1, 0.5)
      Y = rbinom(m, 1, (1-alpha)*3/4 + alpha/4) 
      
      p_M_alpha = discrete.DRPT(X, Y,  r = r)
      
      # Store whether the p-value is below 0.05 (for power calculation)
      pow_alpha[M] = (p_M_alpha < 0.05)
    }
    
    # Store mean and standard deviation of power for this alpha
    pow_vector[i] = mean(pow_alpha)
  }
  
  # Save the results for this `m` value into the list
  all_results[[paste0("m_", m)]] = data.frame(Alpha = alpha_values, Power = pow_vector)
  
  # Plot
  lines(alpha_values, pow_vector, type = "b", col = ind, lwd = 2, pch = 15 + ind)
  ind = ind + 1
}

# Save all data to a CSV file
for (name in names(all_results)) {
  write.csv(all_results[[name]], paste0("experiments/results/", name, "_binary_data.csv"), row.names = FALSE)
}


# Add a horizontal line at 0.05 (significance level) + legend
legend("topleft", inset = 0.05,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($\tau = 1$)'), TeX(r'($\tau = 0.5$)'), 
                  TeX(r'($\tau = 0.2$)'), TeX(r'($\tau = 0.05$)')),
       col = c(1, 2, 3, 4),
       pch = c(15, 16, 17, 18))
abline(h = 0.05, col = "red", lty = 2)
dev.off()