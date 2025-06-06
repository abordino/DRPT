source(file = "./DRPT.R")
source(file = "./estimate_marginal_ratio3.R")

DRPTcond2sample = function(d1, d2, S, H, kernel, prop = 0.5, 
                           est.method = "LL", alpha = 0.05){
  
  #------------------------------------------------------
  # split the data according to the proportion given by prop
  #------------------------------------------------------

  
  n1 = length(d1$y)
  n2 = length(d2$y)
  d = if (is.null(dim(d1$x))) 1 else dim(d1$x)[2]
  
  x1 = d1$x; y1 = as.matrix(d1$y)
  x2 = d2$x; y2 = as.matrix(d2$y)
  
  n1 = length(y1); n2 = length(y2)
  n12 = ceiling(n1 * prop); n22 = ceiling(n2 * prop)
  n11 = n1 - n12; n21 = n2 - n22
  x11 = x1[1:n11, , drop=F]; x12 = x1[-(1:n11), , drop=F]
  y11 = y1[1:n11, , drop=F]; y12 = y1[-(1:n11), , drop=F]
  x21 = x2[1:n21, , drop=F]; x22 = x2[-(1:n21), , drop=F]
  y21 = y2[1:n21, , drop=F]; y22 = y2[-(1:n21), , drop=F]
    
  #------------------------------------------------------
  # Estimate marginal density ratio based on first half of the sample
  #------------------------------------------------------
  r.hat0 = estimate_marginal_ratio3(x11, x21, y11, y21, est.method)
  
  x_values = data.frame(rbind(x12,x22))
  colnames(x_values) = c("V1", "V2", "V3", "V4", "V5", "V6")
  r.values = r.hat0(as.matrix(x_values))
  
  make_lookup_vector = function(x1, y1, r1, x2, y2, r2) {
    
    keys1 = apply(cbind(x1, y1), 1, paste, collapse = "_")
    keys2 = apply(cbind(x2, y2), 1, paste, collapse = "_")
    
    lookup_vec = setNames(c(r1, r2), c(keys1, keys2))
    return(lookup_vec)
  }
  
  lookup_vec = make_lookup_vector(x12, y12, r.values[1:n12],
                                  x22, y22, r.values[(n12+1):(n12+n22)])
  
  r.hat = function(...) {
    z = cbind(...)
    if (!is.matrix(z)) z = matrix(z, ncol = length(z))
    keys = apply(z, 1, paste, collapse = "_")
    return(lookup_vec[keys])
  }
    
  #------------------------------------------------------
  # Run the DRPT on the second half of the sample, using marg_ratio as a parameter
  #------------------------------------------------------
  test0 = cbind(x12,y12); test1 = cbind(x22,y22)
  p.val = DRPT(test0, test1, S, H, r.hat, kernel)
  
  return(p.val)
  
}
