rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
library(Rcpp)
library(future.apply)
library(rootSolve)
library(kernlab)
library(purrr)

sourceCpp("DRPT/shiftedMMD.cpp")
sourceCpp("DRPT/star_sampler.cpp")
source("SPT.R")

SPT_R = function(X, Y, S, H, r, kernel) {
  
  # Calculate sizes of X and Y once
  X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  
  # Generate the sample data using star_sampler
  data = star_sampler_R(X, Y, S, H, r)
  
  # Compute the T_hat_0 value for the original data
  T_hat_0 = shifted.MMD_R(X = data[[1]][1:n,], Y = data[[1]][(n + 1):(n + m),], r, kernel)
  
  # Initialize sum_indicator
  sum_indicator = 0
  
  # Function to calculate shifted MMD for a given sample index
  compute_shifted_MMD = function(b) {
    X_b = data[[1 + b]][1:n,]
    Y_b = data[[1 + b]][(n + 1):(n + m),]
    current_T_hat = shifted.MMD_R(X = X_b, Y = Y_b, r, kernel)
    return(current_T_hat)
  }
  
  # Plan for parallel execution
  plan(multisession)  # Use multisession for parallel execution on local machine
  
  # Use future_sapply for parallel processing
  T_hats = future_sapply(seq_len(H), compute_shifted_MMD)
  
  # Stop the plan to free resources (optional)
  plan(sequential)  # Switch back to sequential processing
  
  # Calculate sum_indicator
  sum_indicator = sum(T_hats >= T_hat_0)
  
  # Compute the final p_hat value
  p_hat = (1 + sum_indicator) / (H + 1)
  return(p_hat)
}

shifted.MMD_R = function(X, Y, r, kernel) {
  X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  tau = n / m
  k = kernel
  
  # Compute lambda.star using root finding
  Z = rbind(X,Y)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k,,drop = FALSE]
      # sum = sum + 1 / (n + m * l * r(Zk))
      sum = sum + 1 / (n + m * l * apply(Zk, 1, function(row) do.call(r, as.list(row))) )
    }
    return(sum - 1)
  }
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
  
  # First term
  first_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in (1:n)[-i]) {
        Xi = X[i,,drop = FALSE]; Xj = X[j,,drop = FALSE]
        r_i =  apply(Xi, 1, function(row) do.call(r, as.list(row)))
        r_j =  apply(Xj, 1, function(row) do.call(r, as.list(row)))
        den = (tau + lambda.star * r_i) * (tau + lambda.star * r_j)
        num = k(Xi, Xj) * r_i * r_j * lambda.star^2
        sum = sum + num / den
      }
    }
    return(sum / (n * (n - 1)))
  }
  
  # Second term
  second_term = function() {
    sum = 0
    for (i in 1:m) {
      for (j in (1:m)[-i]) {
        Yi = Y[i,,drop = FALSE]; Yj = Y[j,,drop = FALSE]
        r_i =  apply(Yi, 1, function(row) do.call(r, as.list(row)))
        r_j =  apply(Yj, 1, function(row) do.call(r, as.list(row)))
        den = (tau + lambda.star * r_i) * (tau + lambda.star * r_j)
        num = k(Yi, Yj)
        sum = sum + num / den
      }
    }
    return(sum / (m * (m - 1)))
  }
  
  # Mixed term
  mixed_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in 1:m) {
        Xi = X[i,,drop = FALSE]; Yj = Y[j,,drop = FALSE]
        r_i = apply(Xi, 1, function(row) do.call(r, as.list(row)))
        r_j =  apply(Yj, 1, function(row) do.call(r, as.list(row)))
        den = (tau + lambda.star * r_i) * (tau + lambda.star * r_j)
        num = k(Xi, Yj) * r_i * lambda.star
        sum = sum + num / den
      }
    }
    return(sum / (n * m))
  }
  
  # Compute and return the shifted MMD
  return(first_term() + second_term() - 2 * mixed_term())
}

star_sampler_R = function(X, Y, S, H, r) {
  X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  h = min(n, m)
  
  Z = rbind(as.matrix(X), as.matrix(Y))
  Z_original = Z
  
  # Define a helper function for swapping based on odds ratio
  swap_based_on_odds_ratio = function(Z, ind_n, ind_m, r) {
    
    odds_ratios = apply(as.matrix(Z[ind_n,]), 1, function(row) do.call(r, as.list(row))) / 
      apply(as.matrix(Z[ind_m,]), 1, function(row) do.call(r, as.list(row)))
    swap_decisions = rbinom(length(ind_n), 1, odds_ratios / (1 + odds_ratios))
    # if (any(is.na(swap_decisions))) break
    
    # Perform swaps where rbinom results in TRUE
    swap_indices = which(swap_decisions == 1)
    if (length(swap_indices) > 0) {
      Z_swap = Z[ind_n[swap_indices],]
      Z[ind_n[swap_indices],] = Z[ind_m[swap_indices],]
      Z[ind_m[swap_indices],] = Z_swap
    }
    
    return(Z)
  }
  
  # Main sampling loop
  for (j in 1:S) {
    ind_n = sample(1:n, h)
    ind_m = sample((n+1):(n+m), h)  # Avoid adding n in the loop
    Z = swap_based_on_odds_ratio(Z, ind_n, ind_m, r)
  }
  
  Z_star = Z
  exch_Z = list(Z_original)  # Store Z_original in the first position
  
  # Bootstrap procedure
  for (b in 1:H) {
    for (j in 1:S) {
      ind_n = sample(1:n, h)
      ind_m = sample((n+1):(n+m), h)
      Z = swap_based_on_odds_ratio(Z, ind_n, ind_m, r)
    }
    
    exch_Z[[b + 1]] = Z  # Store the updated Z in the list
    Z = Z_star  # Reset Z to Z_star for the next iteration
  }
  
  return(exch_Z)
}



#----------------------------------------------------------------------------
# 1D
#----------------------------------------------------------------------------
n = 250
m = 250
d = 1
H = 19

r = function(x) {
  return(exp(-4 * x^2))
}

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}


for(seed in 1:10){
  
  #create data
  set.seed(seed)
  alpha = 0.12
  X = as.matrix(rnorm(n, 0, 1))
  Y = as.matrix(ifelse(rbinom(m, 1, 1/(1 + alpha)) == 1, 
             rnorm(m, 0, 1/3), 
             rexp(m, 1)))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT(X, Y, S = 100, H, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("SPT in C = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_R(X, Y, S = 100, H, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("SPT in R = ", time.taken, "p_val = ", p))
  
  print("-------------------------------")
  
}


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

print("--------------2D now--------------")

for(seed in 1:10){
  #generate data
  set.seed(seed)
  X = cbind(runif(n, 0, 1), runif(n, 0, 1))
  Y = matrix(0, nrow = m, ncol = 2)
  
  alpha = 0
  binary_vec = rbinom(m, 1, 1 / (1 + alpha))
  
  Y[binary_vec == 1, ] = cbind(invCDF(runif(sum(binary_vec))), invCDF(runif(sum(binary_vec))))
  Y[binary_vec == 0, ] = cbind(rbeta(sum(!binary_vec), 0.5, 0.5), rbeta(sum(!binary_vec), 0.5, 0.5))
  
  # run tests, and should give the same p-value??
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_R(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("SPT in C = ", time.taken, "p_val = ", p))
  
  start.time = Sys.time()
  set.seed(seed)
  p = SPT_R(X, Y, S = 100, H = 99, r, gaussian.kernel)
  end.time = Sys.time()
  time.taken = end.time - start.time
  time.taken
  print(paste("SPT in R = ", time.taken, "p_val = ", p))
  
  print("-------------------------------")
}


