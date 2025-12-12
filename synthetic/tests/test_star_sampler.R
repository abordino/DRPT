rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(Rcpp)

sourceCpp("DRPT/star_sampler.cpp")

star_samplerMOINE = function(X, Y, S, H, r) {
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

star_samplerMOINEjust1D = function(X, Y, S, H, r) {
  # X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  h = min(n, m)
  
  # Z = rbind(as.matrix(X), as.matrix(Y))
  Z = rbind(X,Y)
  Z_original = Z
  
  # Define a helper function for swapping based on odds ratio
  swap_based_on_odds_ratio = function(Z, ind_n, ind_m, r) {
    odds_ratios = r(Z[ind_n,]) / r(Z[ind_m,])
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
n = 50
m = 50
h = min(n,m)
H = 1

r = function(x){
  return(1/2+x)
}

invCDF = function(x){
  return((-1+sqrt(1+8*x))/2)
}

X = as.matrix(runif(n, 0, 1), ncol = 1)
Y = runif(m, 0, 1)
Y = as.matrix(invCDF(X), ncol = 1)

set.seed(98)
start.time <- Sys.time()
data1 = star_sampler_C(X, Y, S = 100, H, r)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

set.seed(98)
start.time <- Sys.time()
data2 = star_samplerMOINE(X, Y, S = 100, H, r)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

for(i in 1:length(data1)){
  print(data1[[i]]-data2[[i]])
}

set.seed(98)
start.time <- Sys.time()
data = star_samplerMOINEjust1D(X, Y, S = 100, H = 2, r)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


par(mfrow = c(2,2))
for(h in 2:2){
  hist(data1[[1]][1:n], probability = TRUE, col = "red", xlab = " ")
  hist(data1[[1]][(n+1):(n+m)], probability = TRUE, col = "lightblue", xlab = "")
  hist(data1[[h]][1:n], probability = TRUE, col = "red", xlab = " ")
  hist(data1[[h]][(n+1):(n+m)], probability = TRUE, col = "lightblue", xlab = "")
}


#----------------------------------------------------------------------------
# 2D
#----------------------------------------------------------------------------
n = 300000
m = 300000

# Shift function
r = function(x,y) {
  return(4*x*y)
}

invCDF = function(x){
  return(sqrt(x))
}

X = cbind(runif(n, 0, 1), runif(n, 0, 1))

# Pre-allocate the matrix Y with two columns and m rows
Y <- matrix(0, nrow = m, ncol = 2)

# Generate the random binary vector
alpha = 0
binary_vec <- rbinom(m, 1, 1 / (1 + alpha))

# Fill Y based on binary_vec
Y[binary_vec == 1, ] <- cbind(invCDF(runif(sum(binary_vec))), invCDF(runif(sum(binary_vec))))
Y[binary_vec == 0, ] <- cbind(rbeta(sum(!binary_vec), 0.5, 0.5), rbeta(sum(!binary_vec), 0.5, 0.5))


start.time <- Sys.time()
data = star_sampler_C(X, Y, S = 100, H = 2, r)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
data = star_samplerMOINE(X, Y, S = 100, H = 2, r)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

### sanity check
par(mfrow = c(2,2))
for(h in 2:3){
  hist(data[[1]][1:n, 2], probability = TRUE, col = "red", xlab = " ")
  hist(data[[1]][(n+1):(n+m), 2], probability = TRUE, col = "lightblue", xlab = "")
  hist(data[[h]][1:n, 2], probability = TRUE, col = "red", xlab = " ")
  hist(data[[h]][(n+1):(n+m), 2], probability = TRUE, col = "lightblue", xlab = "")
}

### sanity check
par(mfrow = c(2,2))
for(h in 2:3){
  hist(data[[1]][1:n, 1], probability = TRUE, col = "red", xlab = " ")
  hist(data[[1]][(n+1):(n+m), 1], probability = TRUE, col = "lightblue", xlab = "")
  hist(data[[h]][1:n, 1], probability = TRUE, col = "red", xlab = " ")
  hist(data[[h]][(n+1):(n+m), 1], probability = TRUE, col = "lightblue", xlab = "")
}