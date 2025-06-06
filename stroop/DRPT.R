library(Rcpp)
library(RcppParallel)
library(future.apply)
library(furrr)

sourceCpp("DRPT/shiftedMMD.cpp")
sourceCpp("DRPT/star_sampler.cpp")

DRPT = function(X, Y, S, H, r, kernel) {
  
  X = as.matrix(X); Y = as.matrix(Y)
  n = nrow(X); m = nrow(Y)
  
  # Generate exchangeable dataset using star_sampler
  data = star_sampler_C(X, Y, S, H, r)
  
  # Compute T_hat_0
  T_hat_0 = compute_mmd_C(X = as.matrix(data[[1]][1:n,]), 
                          Y = as.matrix(data[[1]][(n + 1):(n + m),]), 
                          r, kernel)
  
  # Compute T_hat_b
  sum_indicator = 0
  
  compute_shifted_MMD = function(b) {
    X_b = as.matrix(data[[1 + b]][1:n,])
    Y_b = as.matrix(data[[1 + b]][(n + 1):(n + m),])
    current_T_hat = compute_mmd_C(X = X_b, Y = Y_b, r, kernel)
    return(current_T_hat)
  }
  
  plan(multicore) 
  T_hats = future_sapply(seq_len(H), compute_shifted_MMD, future.seed = TRUE)
  
  plan(sequential)  
  sum_indicator = sum(T_hats >= T_hat_0)
  
  p_hat = (1 + sum_indicator) / (H + 1)
  return(p_hat)
}
