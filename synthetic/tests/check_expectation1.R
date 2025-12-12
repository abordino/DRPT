rm(list = ls())  # Clear environment
gc()             # Free memory
.rs.restartR()   # Restart R session (for RStudio users)

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(rootSolve)
library(Rcpp)

source("star_sampler.R")
source("shiftedMMD.R")

sourceCpp("DRPT/shiftedMMD.cpp")


#----------------------------------------------------------------------------
# 1D
#----------------------------------------------------------------------------
n = 300
m = 300
tau = n/m

r = function(x){
  return(1/2+x)
}

sigma = 1

gaussian.kernel = function(x, y){
  return(sigma * exp(-sigma^2*sum(((x - y) ^ 2))))
}

invCDF = function(x){
  return((-1+sqrt(1+8*x))/2)
}

compute_lambda = function(X, Y, r){
  # Compute lambda.star using root finding
  Z = rbind(X,Y)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k]
      sum = sum + 1 / (n + m * l * r(Zk) )
    }
    return(sum - 1)
  }
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^5)[1]
  return(lambda.star)
}

compute.star = function(X, Y, r, k, lambda){
  
  # First term
  first_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in 1:n) {
        Xi = X[i]; Xj = X[j]
        r_i = r(Xi)
        r_j =  r(Xj)
        den = (tau + lambda * r_i) * (tau + lambda * r_j)
        num = k(Xi, Xj) * r_i * r_j * lambda^2
        sum = sum + num / den
      }
    }
    return(sum / n^2)
  }
  
  # Second term
  second_term = function() {
    sum = 0
    for (i in 1:m) {
      for (j in 1:m) {
        Yi = Y[i]; Yj = Y[j]
        r_i =  r(Yi)
        r_j =  r(Yj)
        den = (tau + lambda * r_i) * (tau + lambda * r_j)
        num = k(Yi, Yj)
        sum = sum + num / den
      }
    }
    return(sum / m^2)
  }
  
  # Mixed term
  mixed_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in 1:m) {
        Xi = X[i]; Yj = Y[j]
        r_i = r(Xi)
        r_j = r(Yj)
        den = (tau + lambda * r_i) * (tau + lambda * r_j)
        num = k(Xi, Yj) * r_i * lambda
        sum = sum + num / den
      }
    }
    return(sum / (n * m))
  }
  
  # compute the term S
  S = function() {
    sum = 0
    for (i in 1:n) {
        r_i = r(X[i])
        den = tau + lambda * r_i
        num = lambda * r_i
        sum = sum + num / den
    }
    return(sum / n)
  }
  
  T.squared = first_term() + second_term() - 2*mixed_term()
  
  I. = (S())^2 * sigma/n
  
  II. = -0.5 * (first_term() + second_term()) / n

  return((S() - 0.5/n ) * T.squared - I. - II.)
}

compute.S = function(X, Y, r, lambda){
  
  # compute the term S
  S = function() {
    sum = 0
    for (i in 1:n) {
      r_i = r(X[i])
      den = tau + lambda * r_i
      num = lambda * r_i
      sum = sum + num / den
    }
    return(sum / n)
  }
  
  return(S())
}

compute.A = function(X, Y, r, lambda){
  
  # First term
  first_term = function() {
    sum = 0
    for (i in 1:n) {
        Xi = X[i]
        r_i = r(Xi)
        den = (tau + lambda * r_i)^2
        num = r_i^2 * lambda^2
        sum = sum + num / den
      }
    return(sum)
  }
  
  # Second term
  second_term = function() {
    sum = 0
      for (j in 1:m) {
        Yj = Y[j]
        r_j =  r(Yj)
        den = (tau + lambda * r_j)^2
        sum = sum + 1 / den
      }
    return(sum)
  }
  
  # compute the term S
  S = function() {
    sum = 0
    for (i in 1:n) {
      r_i = r(X[i])
      den = tau + lambda * r_i
      num = lambda * r_i
      sum = sum + num / den
    }
    return(sum / n)
  }
  
  return(S() * (sigma*S()/n - sigma*first_term()/(n^2) - sigma*second_term()/(m^2)))
}

compute.B = function(X, Y, r, lambda, k){
  
  # First term
  first_term = function() {
    sum = 0
    for (i in 1:n) {
      for (j in (1:n)[-i]) {
        Xi = X[i]; Xj = X[j]
        r_i = r(Xi); r_j =  r(Xj)
        den = (tau + lambda * r_i) * (tau + lambda * r_j)
        num = k(Xi, Xj) * r_i * r_j * lambda^2
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
        Yi = Y[i]; Yj = Y[j]
        r_i =  r(Yi); r_j =  r(Yj)
        den = (tau + lambda * r_i) * (tau + lambda * r_j)
        num = k(Yi, Yj)
        sum = sum + num / den
      }
    }
    return(sum / (m * (m - 1)))
  }
  
  # compute the term S
  S = function() {
    sum = 0
    for (i in 1:n) {
      r_i = r(X[i])
      den = tau + lambda * r_i
      num = lambda * r_i
      sum = sum + num / den
    }
    return(sum / n)
  }
  
  return( (2*S()-1) * first_term() + (2*tau*S()-1) * second_term() )
}


#### compute A
#----------------------------------------------------------------------------
n = 3000
m = 1000
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 4000
res = 0

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)
  
  lambda.star = compute_lambda(X, Y, r)
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  res = res + compute.A(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
                                  r = r, lambda = lambda.star)/B
}

print(paste("A.average = ", res))

#### compute B
#----------------------------------------------------------------------------
n = 3000
m = 1000
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 4000
res = 0

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)
  
  lambda.star = compute_lambda(X, Y, r)
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  res = res + compute.B(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
                        r = r, lambda = lambda.star, k = gaussian.kernel)/B
}

print(paste("B.average = ", res*B/b))

#### compute correlation
#----------------------------------------------------------------------------
n = 600
m = 200
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 1000
S.vec = c()
T.squared.vec = c()

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)
  
  lambda.star = compute_lambda(X, Y, r)
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  
  S.vec = c(S.vec, compute.S(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
                        r = r, lambda = lambda.star))
  T.squared.vec = c(T.squared.vec, 
                    shifted.MMD(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
                                r = r, k = gaussian.kernel))
}

plot(T.squared.vec, S.vec)
print(paste("Corr(S, T^2) = ", cor(T.squared.vec, S.vec)))




#### S_n
#----------------------------------------------------------------------------
n = 3000
m = 1000
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 4000
res.perm = 0
res.id = 0
alpha = 0.8

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  
  mixture_component = runif(m)
  invCDF_samples = invCDF(runif(m, 0, 1))
  beta_samples = rbeta(m, 0.05, 1)
  Y = as.matrix(ifelse(mixture_component < alpha, beta_samples, invCDF_samples), 
                ncol = 1)
  
  res.id = res.id + compute.S(X = X, Y = Y, r = r, lambda = lambda.star)/B
  
  lambda.star = compute_lambda(X, Y, r)
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  res.perm = res.perm + compute.S(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
                           r = r, lambda = lambda.star)/B
}

print(paste("S.id = ", res.id))
print(paste("S.perm = ", res.perm))
print(paste("1/1+tau = ", 1/(1+tau)))


#### inspect S for permuted data (similar to the null)
n = 100000
m = 200000
tau = n/m

r = function(x){
  return(as.numeric(0*x + 1))
}

invCDF = function(x){
  return((-1+sqrt(1+8*x))/2)
}

X = as.matrix(runif(n, 0, 1), ncol = 1)
Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)

lambda.star = compute_lambda(X, Y, r)

permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]

compute.S(X = permuted.data[1:n], Y = permuted.data[(n+1):(n+m)],
          r = r, lambda = lambda.star)





#### compute exp T^2
#----------------------------------------------------------------------------
n = 500
m = 300
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 10000
S.vec = c()
T.squared.vec_01 = c()

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)
  
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  
  T.squared.vec_01 = c(T.squared.vec_01, 
                       compute_mmd_C(X = as.matrix(permuted.data[1:n]), 
                                     Y = as.matrix(permuted.data[(n+1):(n+m)]),
                                     r = r, k = gaussian.kernel))
}

print(paste("E(T^2) = ", mean(T.squared.vec_01)))

#### compute exp T^2
#----------------------------------------------------------------------------
n = 200
m = 200
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 10000
S.vec = c()
T.squared.vec_02 = c()

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(invCDF(runif(m, 0, 1)), ncol = 1)
  
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  
  T.squared.vec_02 = c(T.squared.vec_02, 
                       compute_mmd_C(X = as.matrix(permuted.data[1:n]), 
                                     Y = as.matrix(permuted.data[(n+1):(n+m)]),
                                     r = r, k = gaussian.kernel))
}

print(paste("E(T^2) = ", mean(T.squared.vec_02)))



#### compute exp T^2
#----------------------------------------------------------------------------
n = 500
m = 300
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 10000
S.vec = c()
T.squared.vec_11 = c()

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(runif(m, 0, 1), ncol = 1)
  
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  
  T.squared.vec_11 = c(T.squared.vec_11, 
                       compute_mmd_C(X = as.matrix(permuted.data[1:n]), 
                                     Y = as.matrix(permuted.data[(n+1):(n+m)]),
                                     r = r, k = gaussian.kernel))
}

print(paste("E(T^2) = ", mean(T.squared.vec_11)))

#### compute exp T^2
#----------------------------------------------------------------------------
n = 200
m = 200
tau = n/m

r = function(x){
  return(1/2 + x)
}

B = 10000
S.vec = c()
T.squared.vec_12 = c()

for (b in 1:B){
  print(b)
  X = as.matrix(runif(n, 0, 1), ncol = 1)
  Y = as.matrix(runif(m, 0, 1), ncol = 1)
  
  permuted.data = star_sampler(X,Y,S = 200, H = 1, r)[[2]]
  
  T.squared.vec_12 = c(T.squared.vec_12, 
                       compute_mmd_C(X = as.matrix(permuted.data[1:n]), 
                                     Y = as.matrix(permuted.data[(n+1):(n+m)]),
                                     r = r, k = gaussian.kernel))
}

print(paste("E(T^2) = ", mean(T.squared.vec_12)))


