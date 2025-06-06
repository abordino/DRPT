library(CVST)

estimate_marginal_ratio3 = function(x1, x2, y1, y2, est.method="LL", seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n1 = length(y1); n2 = length(y2)
  label.fit = factor(c(rep(0,n1), rep(1,n2)))
  
  if (est.method == "LL"){
    x.fit = rbind(x1,x2)
    fit.marginal = glm(label.fit~., data=as.data.frame(x.fit), family=binomial())
    
    marginal_ratio = function(new_data){
      prob.marginal = predict(fit.marginal, newdata=as.data.frame(new_data), type="response")
      ratio = prob.marginal/(1-prob.marginal)*n1/n2
      return(as.numeric(ratio))
    }
    
    return(marginal_ratio)
    
  } else if (est.method == "KLR"){
    
    klrlearner = constructKlogRegLearner()
    params = list(kernel='rbfdot', sigma=0.005, lambda=0.0005, tol=1e-6, maxiter=500)
    
    x.fit = rbind(x1, x2)
    data.fit = constructData(x.fit, label.fit)
    fit.marginal = klrlearner$learn(data.fit, params)
    
    marginal_ratio = function(newdata){
      K = kernelMult(fit.marginal$kernel, newdata, 
                     fit.marginal$data, fit.marginal$alpha)
      pi = 1 / (1 + exp(-as.vector(K)))
      ratio = pi/(1-pi)*n1/n2
      
      return(as.numeric(ratio))
    }
    
    return(marginal_ratio)
    
  }
}
