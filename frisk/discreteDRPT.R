library(BiasedUrn)

source(file = "T_hat_discrete.R")

discrete.DRTP = function(X, Y, r, H, type = "D", lambda = 1) {
  lev = sort(union(X,Y))
  NX = table(factor(X, levels = lev)); NY = table(factor(Y, levels = lev))
  tot = NX + NY  
  n = length(X)
  m = length(Y)
  r = r/r[1]
  
  # Compute T_hat_0 
  T_hat_0 = T_hat.discrete(NY, r, n, m, tot, type, lambda)
  
  # Compute T_hat_b
  sum_indicator = 0
  
  for (b in 1:H) {
    N_sigma = rMFNCHypergeo(1, tot, m, r)
    T_hat_b = T_hat.discrete(N_sigma, r, n, m, tot, type, lambda)
    
    sum_indicator = sum_indicator + (T_hat_b >= T_hat_0)
  }
  
  p_hat = (1 + sum_indicator) / (H + 1)
  return(p_hat)
}
