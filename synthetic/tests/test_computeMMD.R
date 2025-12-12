rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
library(Rcpp)
library(rootSolve)

sourceCpp("DRPT/shiftedMMD.cpp")

shifted.MMD_MOINE = function(X, Y, r, kernel) {
  # X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  tau = n / m
  k = kernel
  
  # Compute lambda.star using root finding
  Z = rbind(X,Y)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k, ,drop = FALSE]
      sum = sum + 1 / (n + m * l * apply(Zk, 1, function(row) do.call(r, as.list(row))) )
    }
    return(sum - 1)
  }
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)
  
  # First term
  first_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in (1:n)[-i]) {
        Xi = X[i,,drop = FALSE]; Xj = X[j,,drop = FALSE]
        r_i = apply(Xi, 1, function(row) do.call(r, as.list(row)))
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
        r_j = apply(Yj, 1, function(row) do.call(r, as.list(row)))
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


#----------------------------------------------------------------------------
# 1D - uniform
#----------------------------------------------------------------------------
n = 300
m = 300

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


start.time <- Sys.time()
mmd1 = compute_mmd_C(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
mmd2 = shifted.MMD_MOINE(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

mmd1
mmd2

#----------------------------------------------------------------------------
# 1D - gaussian
#----------------------------------------------------------------------------
n = 200
m = 200

r = function(x) {
  return(exp(-4 * x^2))
}

gaussian.kernel = function(x, y, lambda = 1){
  return(exp(-sum(((x - y) ^ 2) / (2* lambda ^ 2))))
}

  
alpha = 0.0
X = as.matrix(rnorm(n, 0, 1))
Y = as.matrix(ifelse(rbinom(m, 1, 1/(1 + alpha)) == 1, 
                       rnorm(m, 0, 1/3), 
                       rexp(m, 1)))

start.time <- Sys.time()
mmd1 = compute_mmd_C(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
mmd2 = shifted.MMD_MOINE(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

mmd1
mmd2


#----------------------------------------------------------------------------
# 2D
#----------------------------------------------------------------------------
n = 300
m = 300

# Shift function
r = function(x,y) {
  return(4*x*y)
}

# Gaussian kernel 
gaussian.kernel <- function(x, y = NULL, h = 1) {
  if (is.null(y)) {
    res <- (exp(-0.5 * (x/h)^2) / (h * sqrt(2 * pi)))
  } else {
    dist <- sum((x - y)^2)
    res <- exp(-dist / (2 * h^2))
  }
  return(res)
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
mmd1 = compute_mmd_C(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
mmd2 = shifted.MMD_MOINE(X, Y, r, gaussian.kernel)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

mmd1
mmd2
