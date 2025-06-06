T_hat.discrete = function(N, r, n, m, tot, type = "D", lambda = 1) {
  K = length(N)
  
  if (type == "V"){
    denominator = (n / m) + lambda * r[1:K]
    f_j = (tot[1:K] - N[1:K]) / n
    g_j = N[1:K] / m
    
    term1 = (lambda * r[1:K] * f_j) / denominator
    term2 = g_j / denominator
    
    return(sum((term1 - term2)^2))
  }
  
  else if (type == "U"){
    denominator = (n / m) + lambda * r[1:K]
    f_j = (tot[1:K] - N[1:K]) / n
    g_j = N[1:K] / m
    
    term1 = (lambda * r[1:K] * f_j) / denominator
    term2 = g_j / denominator
    
    V = sum((term1 - term2)^2)
    
    termX = sum(((lambda * r[1:K]) / denominator)^2 * ((tot[1:K] - N[1:K]) / n^2))
    termY = sum((1 / denominator)^2 * (N[1:K] / m^2))
    
    return(V - termX - termY)
  }
  
  else if (type == "D"){
    term1 = as.numeric(N[2:K]) * as.numeric((tot[1] - N[1])) / as.numeric(sqrt(r[2:K]))
    term2 = as.numeric(sqrt(r[2:K])) * as.numeric(tot[2:K] - N[2:K]) * as.numeric(N[1])
    
    return(sum(abs(term1 - term2)) / (n * m))
  }
  
}
