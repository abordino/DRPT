rm(list = ls())  # Clear environment
gc()             # Free memory
# .rs.restartR()   # Restart R session (for RStudio users)

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
library(Rcpp)
library(RcppParallel)
library(future.apply)
library(furrr)

sourceCpp("DRPT/shiftedMMD.cpp")
sourceCpp("DRPT/star_sampler.cpp")
sourceCpp("DRPT/SPT.cpp")
# sourceCpp("DRPT/SPT_par.cpp")


SPT_MOINE = function(X, Y, S, H, r, kernel) {
  # Ensure matrices and dimensions
  X = as.matrix(X); Y = as.matrix(Y)
  n = nrow(X); m = nrow(Y)
  
  # Generate exchangeable datasets using star_sampler
  data = star_sampler_C(X, Y, S, H, r)
  
  # Compute T_hat_0 for the original dataset
  T_hat_0 = compute_mmd_C(X = as.matrix(data[[1]][1:n,]), 
                           Y = as.matrix(data[[1]][(n + 1):(n + m),]), 
                           r, kernel)
  
  # Vectorized calculation of shifted MMD for bootstrap samples
  T_hats = sapply(1:H, function(b) {
    X_b = as.matrix(data[[1 + b]][1:n, ])
    Y_b = as.matrix(data[[1 + b]][(n + 1):(n + m), ])
    compute_mmd_C(X = X_b, Y = Y_b, r, kernel)
  })
  
  # Calculate sum_indicator and p_hat
  sum_indicator = sum(T_hats >= T_hat_0)
  p_hat = (1 + sum_indicator) / (H + 1)
  
  return(p_hat)
}

SPT_MOINE_par = function(X, Y, S, H, r, kernel) {
  # Ensure matrices and dimensions
  X = as.matrix(X); Y = as.matrix(Y)
  n = nrow(X); m = nrow(Y)
  
  # Generate exchangeable datasets using star_sampler
  data = star_sampler_C(X, Y, S, H, r)
  
  # Compute T_hat_0 for the original dataset
  T_hat_0 = compute_mmd_C(X = as.matrix(data[[1]][1:n,]), 
                           Y = as.matrix(data[[1]][(n + 1):(n + m),]), 
                           r, kernel)
  
  # Use all available cores except 1
  plan(multicore, workers = min(H, parallel::detectCores() - 1))
  
  # Parallel computation of shifted MMD using future_sapply
  T_hats = future_sapply(1:H, function(b) {
    X_b = as.matrix(data[[1 + b]][1:n, ])
    Y_b = as.matrix(data[[1 + b]][(n + 1):(n + m), ])
    compute_mmd_C(X = X_b, Y = Y_b, r, kernel)
  }, future.seed = TRUE)  # Ensure reproducibility
  
  # Reset to sequential to avoid memory leaks
  plan(sequential)
  
  # Calculate p_hat
  sum_indicator = sum(T_hats >= T_hat_0)
  p_hat = (1 + sum_indicator) / (H + 1)
  
  return(p_hat)
}





#----------------------------------------------------------------------------
# 1D
#----------------------------------------------------------------------------
n = 30
m = 30

r = function(x){
  return(1/2+x)
}

gaussian.kernel = function(x, y, lambda = 1){
  return(exp(-sum(((x - y) ^ 2) / (2* lambda ^ 2))))
}

invCDF = function(x){
  return((-1+sqrt(1+8*x))/2)
}

X = as.matrix(runif(n, 0, 1), ncol = 1)
Y = as.matrix(runif(m, 0, 1))
Y = as.matrix(invCDF(X), ncol = 1)

for(seed in 1:10){

  start.time = Sys.time()
  set.seed(seed)
  p = SPT_C(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("Rccp_C = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_MOINE(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("Rccp-Pingpong = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_MOINE_par(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("Rccp+par = ", time.taken, "p_val = ", p))
  
  print("-------------------------------")

}

# start.time = Sys.time()
# set.seed(seed)
# SPT_C_par(X, Y, S = 100, H = 99, r, gaussian.kernel)
# end.time = Sys.time()
# time.taken = end.time - start.time
# time.taken

#----------------------------------------------------------------------------
# 2D
#----------------------------------------------------------------------------
n = 30
m = 30

r = function(x,y){
  return(4*x*y)
}

gaussian.kernel = function(x, y, lambda = 1){
  return(exp(-sum(((x - y) ^ 2) / (2* lambda ^ 2))))
}

# Shift function
r = function(x,y) {
  return(4*x*y)
}

invCDF = function(x){
  return(sqrt(x))
}

X = cbind(runif(n, 0, 1), runif(n, 0, 1))

# Pre-allocate the matrix Y with two columns and m rows
Y = matrix(0, nrow = m, ncol = 2)

# Generate the random binary vector
alpha = 0
binary_vec = rbinom(m, 1, 1 / (1 + alpha))

# Fill Y based on binary_vec
Y[binary_vec == 1, ] = cbind(invCDF(runif(sum(binary_vec))), invCDF(runif(sum(binary_vec))))
Y[binary_vec == 0, ] = cbind(rbeta(sum(!binary_vec), 0.5, 0.5), rbeta(sum(!binary_vec), 0.5, 0.5))

print("--------------2D now--------------")

for(seed in 1:10){

  start.time = Sys.time()
  set.seed(seed)
  p = SPT_C(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("Rccp_C = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_MOINE(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("Rccp-Pingpong = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_MOINE_par(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  print(paste("Rccp+par = ", time.taken, "p_val = ", p))

  print("-------------------------------")
}

# start.time = Sys.time()
# set.seed(seed)
# SPT_C_par(X, Y, S = 100, H = 99, r, gaussian.kernel)
# end.time = Sys.time()
# time.taken = end.time - start.time
# time.taken


