rm(list = ls())  # Clear environment
gc()             # Free memory

# setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
setwd("/storage/stats/strtng/SPT/synthetic")

library(MASS)
library(latex2exp)
library(squash)

source(file = "DRPT.R")

##################################################################
# Set parameters
##################################################################
n = 250
m = 250
d = 1

MC = 200
S = 50

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

invCDF = function(x){
  return(sqrt(x))
}


##################################################################
# Simulation DRPT using LL estimator instead of the true r_star
##################################################################
size_n = numeric(4)
n_index = 1
alpha = 0

for (n_train in c(100, 200, 1000, 2000)){
  X = runif(n_train, 0, 1)
  Y = invCDF(runif(n_train, 0, 1))
  n_x <- length(X)
  n_y <- length(Y)
  
  df_train <- data.frame(x = c(X, Y), label  = factor(c(rep(0, n_x), rep(1, n_y))))
  fit.marginal <- glm(label ~ x, data = df_train, family = binomial())
  
  marginal_ratio <- function(z){
    prob <- predict(fit.marginal,newdata = data.frame(x = z), type = "response")         
    return(prob / (1 - prob))                      
  }
  
  sum_ind = 0
  
  # run the simulation
  for (b in 1:MC) {
    # generate data
    X = runif(n, 0, 1)
    Y = invCDF(runif(m, 0, 1))
    
    # define r.hat
    r.values = marginal_ratio(c(X,Y))
    
    create_lookup_list = function(X, rX, Y, rY) {
      keysX = paste(X)
      keysY = paste(Y)
      
      lookup_list = setNames(as.list(c(rX, rY)), c(keysX, keysY))
      
      return(lookup_list)
    }
    
    lookup_table = create_lookup_list(X, r.values[1:n], Y, tail(r.values, m))
    
    r.hat = function(...) {
      z = cbind(...)
      if (is.matrix(z) || is.data.frame(z)) {
        keys = apply(z, 1, paste, collapse = "_")
        return(sapply(keys, function(key) ifelse(key %in% names(lookup_table), lookup_table[[key]], NA)))
      } else {
        key = paste(z, collapse = "_")
        return(ifelse(key %in% names(lookup_table), lookup_table[[key]], NA))
      }
    }
    
    # Our test: call the DRPT function
    p_M_alpha = DRPT(X, Y, S, H = 99, r.hat, kernel = gaussian.kernel)
    
    if (p_M_alpha < 0.05) {
      sum_ind = sum_ind + 1
    }
  }
  
  size_n[n_index] = sum_ind / MC
  print(sum_ind / MC)
  n_index = n_index + 1
}

######################################################
# Save the result in a .csv
######################################################
data_to_save = data.frame(n = c(100, 200, 1000, 2000), Size = size_n)
write.csv(data_to_save, "experiments/results/r_estimation_LL.csv", row.names = FALSE)
