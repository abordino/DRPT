source("DRPT.R")

rej.sampling.MMD = function(X, Y, S, H, r, kernel, M) {
  
  n = length(X)
  m = length(Y)
  
  W = c()
  for(i in 1:n){
    U = runif(1)
        
    if (U <= r(X[i]) / M) {
      W = c(W, X[i])
    }
  }
  
  if(length(W) > 0){
    p_hat = DRPT(W, Y, S, H, r = function(x) as.numeric(0*x + 1), kernel)
  } else {
    p_hat = 1
  }
 
  return(p_hat)
}
