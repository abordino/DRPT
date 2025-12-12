X = as.matrix(rnorm(100000, 0, 1))
hist(X, probability = TRUE)
M = 1
r = function(x) {
  return(exp(-4 * x^2))
}

W = c()
for(i in 1:100000){
  U = runif(1) # Step 1: Draw a uniform random variable to compare with X_i
  
  # Step 2: Accept x with probability r(x)/M
  if (U <= r(X[i]) / M) {
    W = c(W, X[i])
  }
}

hist(W, probability = TRUE)
curve(dnorm(x, mean = 0, sd = 1/3), 
      col = "red", lwd = 2, add = TRUE)