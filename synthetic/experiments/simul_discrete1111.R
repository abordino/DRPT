rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(latex2exp)

source(file = "discreteDRPT.R")
source(file = "DRPT.R")
source(file = "rej_sampling_MMD.R")

#####################################################################
# set parameters for d = 4, (1,1,1,1), Test statistic D in the Appendix
#####################################################################
n = 250
m = 250
r = c(1, 1, 1, 1)
MC = 1000

##################################################################
# Simulation for varying alpha
##################################################################
pow_vector = numeric(30)
alpha_values = seq(from = 0, to = 0.8, length.out = 30)

for (i in seq_along(alpha_values)) {
  alpha = alpha_values[i]
  gamma = (25*alpha/43) * (1/sqrt(3) + 1/sqrt(10))
  
  pow_alpha = numeric(MC)
  
  for (M in 1:MC) {
    # Sampling for X and Y based on the probabilities
    X = sample(0:3, n, prob = c(1/4, 1/4, 1/4, 1/4), replace = TRUE)
    Y = sample(0:3, m, prob = c(1/4, 1/4, (1+gamma)/4, (1-gamma)/4), 
               replace = TRUE)
    
    # Compute p-value using DRPT
    p_M_alpha = discrete.DRPT(X, Y, r, H = 99) 
    pow_alpha[M] = (p_M_alpha < 0.05)
  }
  
  pow_vector[i] = mean(pow_alpha)
}

#####################################################################
# Save the data to a CSV file
#####################################################################
data_to_save = data.frame(Alpha = alpha_values, Power = pow_vector)
write.csv(data_to_save, "experiments/results/4d_1111_discrete_data.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------------------

#####################################################################
# set parameters for d = 4, (1,1,1,1) - U/V
#####################################################################
n = 250
m = 250
r = c(1, 1, 1, 1)
MC = 3000
r.func = function(x){
  return(r[x+1])
}

##################################################################
# Simulation for varying alpha
##################################################################
pow_vector_U = numeric(15); pow_vector_V = numeric(15)
alpha_values = seq(from = 0, to = 0.4, length.out = 15)

for (i in seq_along(alpha_values)) {
  alpha = alpha_values[i]
  gamma = (25*alpha/43) * (1/sqrt(3) + 1/sqrt(10))
  
  pow_alpha_U = numeric(MC)
  pow_alpha_V = numeric(MC)
  
  for (M in 1:MC) {
    # Sampling for X and Y based on the probabilities
    X = sample(0:3, n, prob = c(1/4, 1/4, 1/4, 1/4), replace = TRUE)
    Y = sample(0:3, m, prob = c(1/4, 1/4, (1+gamma)/4, (1-gamma)/4), 
               replace = TRUE)
    
    # compute hat.lambda for the test statistic
    Z = rbind(as.matrix(X),as.matrix(Y))
    sum_lambda = function(l) {
      sum = 0
      for (k in 1:(n+m)) {
        Zk = Z[k,,drop = FALSE]
        sum = sum + 1 / (n + m * l * r.func(Zk))
      }
      return(sum - 1)
    }
    lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
    
    # Compute p-value using DRPT
    p_M_alpha_U = discrete.DRPT(X, Y, r, H = 99, type = "U", lambda = lambda.star) 
    p_M_alpha_V = discrete.DRPT(X, Y, r, H = 99, type = "V", lambda = lambda.star) 

    pow_alpha_U[M] = (p_M_alpha_U < 0.05)
    pow_alpha_V[M] = (p_M_alpha_V < 0.05)
  }
  
  pow_vector_U[i] = mean(pow_alpha_U)
  pow_vector_V[i] = mean(pow_alpha_V)
}

#####################################################################
# Save the data to a CSV file
#####################################################################
data_to_save = data.frame(Alpha = alpha_values, Power_V = pow_vector_V,  
                           Power_U = pow_vector_U)
write.csv(data_to_save, "experiments/results/4d_1111_discrete_data_UV.csv", row.names = FALSE)

#---------------------------------------------------------------------------------------------------

#####################################################################
# Plot the results
#####################################################################

# data.MMD = read.csv(file = "experiments/results/4d_1111_MMD.csv")
# data.discrete = read.csv(file = "experiments/results/4d_1111_discrete_data.csv")
data.discrete_UV = read.csv(file = "experiments/results/4d_1111_discrete_data_UV.csv")

png("experiments/pictures/4d_1111.png")
par(mfrow = c(1, 1))
plot(data.discrete_UV$Alpha, data.discrete_UV$Power_U, type = "b", col = "darkviolet", lwd = 2,
     pch = 14, xlab = TeX(r'($\eta$)', bold = TRUE),
     ylab = "Power", ylim = c(0,1))
lines(data.discrete_UV$Alpha, data.discrete_UV$Power_V, type = "b", col = "blue", lwd = 2,
     pch = 15)
legend("topleft", inset = +0.01,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("Discrete DRPT - U", "Discrete DRPT - V"),
       col = c("darkviolet", "blue"),
       pch = c(14, 14))
abline(h = 0.05, col = "red", lty = 2)
dev.off()
