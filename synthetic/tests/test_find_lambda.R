rm(list = ls())  # Clear environment
gc()             # Free memory
.rs.restartR()   # Restart R session (for RStudio users)

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(Rcpp)
library(rootSolve)

sourceCpp("DRPT/find_lambda.cpp")

compute_lambda_MOINE = function(X, Y, r){
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  
  # Compute lambda.star using root finding
  Z = rbind(X,Y)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k,,drop = FALSE]
      sum = sum + 1 / (n + m * l * apply(Zk, 1, function(row) do.call(r, as.list(row))) )
    }
    return(sum - 1)
  }
  
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
  return(lambda.star)
}


##_----------------------
############ 1D
##_----------------------
# Example matrices
X = matrix(runif(1000), ncol = 1)
Y = matrix(runif(1000), ncol = 1)

# Example r function
r = function(x,y){
  return(1/2 + x^2)
}


# Find lambda.starstart.time <- Sys.time()
start.time <- Sys.time()
print(compute_lambda_star_C(X, Y, r))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
print(compute_lambda_MOINE(X, Y, r))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken





##_----------------------
############ 2D
##_----------------------
# Example matrices
X = matrix(runif(1000), ncol = 2)
Y = matrix(runif(1000), ncol = 2)

# Example r function
r = function(x,y){
  return(1/2 + x*y)
}

compute_lambda_MOINE = function(X, Y, r){
  n = dim(X)[1]; m = dim(Y)[1]; d = dim(X)[2]
  
  # Compute lambda.star using root finding
  Z = rbind(X,Y)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k,,drop = FALSE]
      sum = sum + 1 / (n + m * l * apply(Zk, 1, function(row) do.call(r, as.list(row))) )
    }
    return(sum - 1)
  }
  
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
  return(lambda.star)
}


# Find lambda.starstart.time <- Sys.time()
start.time <- Sys.time()
print(compute_lambda_star_C(X, Y, r))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
print(compute_lambda_MOINE(X, Y, r))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



