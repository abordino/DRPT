rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(latex2exp)
library(rootSolve)

source(file = "discreteDRPT.R")
source(file = "DRPT.R")
source(file = "rej_sampling_MMD.R")
source(file = "T_hat_discrete.R")

#####################################################################
# setting the parameters for d = 2, varying m setting
#####################################################################

n = 100
MC = 5000
r = c(1, 3)
r.func = function(x){
  return(r[x+1])
}

##################################################################
# Simulation for varying m
##################################################################
all_results = list()

# Preparing the plot
png("experiments/pictures/binary.png")
par(mfrow = c(1, 1))
plot(1, type = "n",  xlim = c(0,0.9), ylim = c(0,1),
     xlab = TeX(r'($\eta$)', bold = TRUE), 
     ylab = "Power")

ind = 1
for (m in c(100, 200, 500, 2000)){

  pow_vector = numeric(30)
  alpha_values = seq(from = 0, to = 0.9, length.out = 30)
  
  for (i in seq_along(alpha_values)) {
    alpha = alpha_values[i]
    pow_alpha = numeric(MC)
    
    for (M in 1:MC) {
      # Sampling for X and Y based on the probabilities
      X = rbinom(n, 1, 0.5)
      Y = rbinom(m, 1, (1-alpha)*3/4 + alpha/4) 
      
      # compute hat.lambda for the test statistic
      Z = rbind(as.matrix(X),as.matrix(Y))
      sum_lambda = function(l) {
        sum = 0
        for (k in 1:(n+m)) {
          Zk = Z[k, ,drop = FALSE]
          sum = sum + 1 / (n + m * l * r.func(Zk))
        }
        return(sum - 1)
      }
      lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
      
      # Compute p-value using DRPT
      p_M_alpha = discrete.DRPT(X, Y, r, H = 99, type = "V", lambda = lambda.star) 
      pow_alpha[M] = (p_M_alpha < 0.05)
    }
  
    pow_vector[i] = mean(pow_alpha)
  }
  
  all_results[[paste0("m_", m)]] = data.frame(Alpha = alpha_values, Power = pow_vector)
  
  # add line to the plot
  lines(alpha_values, pow_vector, type = "b", col = ind, lwd = 2, pch = 15 + ind)
  ind = ind + 1
}

# Save all data to a CSV file
for (name in names(all_results)) {
  write.csv(all_results[[name]], paste0("experiments/results/", name, "_binary_data.csv"), row.names = FALSE)
}

# Conclude the plot
legend("topleft", inset = 0.05,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($\tau = 1$)'), TeX(r'($\tau = 0.5$)'), 
                  TeX(r'($\tau = 0.2$)'), TeX(r'($\tau = 0.05$)')),
       col = c(1, 2, 3, 4),
       pch = c(15, 16, 17, 18))
abline(h = 0.05, col = "red", lty = 2)
dev.off()

#####################################################################
# setting the parameters for d = 2, varying r (APPENDIX)
#####################################################################
n = 500; m = n
MC = 3000

all_results = list()

# Preparing the plot
png("experiments/pictures/binaryRvarying.png")
par(mfrow = c(1, 1))
plot(1, type = "n",  xlim = c(0,0.1), ylim = c(0,1),
     xlab = TeX(r'($\alpha$)', bold = TRUE), 
     ylab = "Power")

ind = 1
for (r2 in c(1/10, 1/2, 1, 2, 10)){
  r = c(1, r2)
  
  pow_vector = numeric(30)
  alpha_values = seq(from = 0, to = 0.1, length.out = 30)
  
  for (i in seq_along(alpha_values)) {
    alpha = alpha_values[i]
    gamma = 2*sqrt(r2)*alpha
    pow_alpha = numeric(MC)
    
    for (M in 1:MC) {
      # Sampling for X and Y based on the probabilities
      X = rbinom(n, 1, 0.5)
      Y = abs(rbinom(m, 1, (r2 + gamma)/(1+r2) ) - 0)
      
      # Compute p-value using DRPT
      p_M_alpha = discrete.DRPT(X, Y, r, H = 99) 
      pow_alpha[M] = (p_M_alpha < 0.05)
    }
    
    pow_vector[i] = mean(pow_alpha)
  }
  
  all_results[[paste0("r_", r2)]] = data.frame(Alpha = alpha_values, Power = pow_vector)
  
  # add line to the plot
  lines(alpha_values, pow_vector, type = "b", col = ind, lwd = 2, pch = 14 + ind)
  ind = ind + 1
  
}

# Save all data to a CSV file
for (name in names(all_results)) {
  write.csv(all_results[[name]], paste0("experiments/results/", name, "_binary_data.csv"), row.names = FALSE)
}

# conclude the plot
legend("topleft", inset = 0.05,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($r = 1/10$)'), TeX(r'($r = 1/2$)'),
                  TeX(r'($r = 1$)'), TeX(r'($r = 2$)'), TeX(r'($r = 10$)')),
       col = c(1, 2, 3, 4, 5),
       pch = c(15, 16, 17, 18, 19))
abline(h = 0.05, col = "red", lty = 2)
dev.off()
