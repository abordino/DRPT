rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
# setwd("/storage/stats/strtng/SPT/synthetic")

library(MASS)
library(latex2exp)

source(file = "DRPT.R")
source(file = "rej_sampling_MMD.R")

#####################################################################
# set parameters
#####################################################################

n = 250
m = 250
d = 1
MC = 350
S = 50
xxx = 15
M = 1

r = function(x) {
  return(exp(-4 * x^2))
}

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

##################################################################
# Simulation for varying alpha
##################################################################
pow_vector = numeric(xxx)
pow_vector.rej = numeric(xxx)
alpha_index = 1

for (alpha in seq(from = 0, to = 0.25, length.out = xxx)){
  print(alpha)
  sum_ind = 0
  sum_ind.rej = 0
  
  # Monte Carlo simulation
  for (b in 1:MC) {
    
    X = as.matrix(rnorm(n, 0, 1))
    
    Y = as.matrix(ifelse(rbinom(m, 1, 1/(1 + alpha)) == 1, 
               rnorm(m, 0, 1/3), 
               rexp(m, 1)))
    
    median.bandwidth = median(as.vector(dist(rbind(X,Y), method = "euclidean", diag = FALSE, upper = FALSE)))
    
    # Our test: call the DRPT function
    p_M_alpha = DRPT(X, Y, S, H = 99, r, kernel = function(x,y) gaussian.kernel(x, y, lambda = median.bandwidth))
    
    if (p_M_alpha < 0.05) {
      sum_ind = sum_ind + 1
    }
    
    # Standard MMD with rejection sampling
    p_M_alpha.rej = rej.sampling.MMD(X, Y, S, H = 99, r,
                                     kernel = function(x,y) gaussian.kernel(x, y, lambda = median.bandwidth), M = M)

    if (p_M_alpha.rej < 0.05) {
      sum_ind.rej = sum_ind.rej + 1
    }
  }
  
  pow_vector[alpha_index] = sum_ind / MC
  pow_vector.rej[alpha_index] = sum_ind.rej / MC
  
  alpha_index = alpha_index + 1
  
}

##################################################################
# Save the data to a CSV file
##################################################################
data_to_save = data.frame(Alpha = seq(from = 0, to = 0.25, length.out = xxx),
                          Power = pow_vector, Power_rej = pow_vector.rej)
write.csv(data_to_save, "experiments/results/simul_shiftedMMD.csv", row.names = FALSE)

##################################################################
# Plotting
##################################################################
data.MMD = read.csv(file = "experiments/results/simul_1D.csv")
data.x2 = read.csv(file = "experiments/results/simul1D_x2.csv")
data.rej = read.csv(file = "experiments/results/simul_rej1D.csv")

png("experiments/pictures/simul_1D.png")
par(mfrow = c(1, 1))
plot(data.MMD$Alpha, data.MMD$Power, type = "b", 
     col = "darkviolet", lwd = 2, pch = 14, 
     xlab = TeX(r'($\eta$)', bold = TRUE), ylab = "Power")
lines(data.x2$Alpha, data.x2$Power, type = "b", 
      col = "orange", lwd = 2, pch = 14)
lines(data.rej$Alpha, data.rej$Power_rej, type = "b", 
      col = "green", lwd = 2, pch = 16)
legend("topleft", inset = 0.05,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("DRPT - exp(-4x^2)", "DRPT - exp(-x^2)",
                  "MMD Rejection - exp(-4x^2)"),
       col = c("darkviolet", "orange", "green"),
       pch = c(14, 14, 16))
abline(h = 0.05, col = "red", lty = 2)
dev.off()