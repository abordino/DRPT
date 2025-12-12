rm(list = ls())
gc()

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")
# setwd("/storage/stats/strtng/SPT/synthetic")

library(MASS)
library(latex2exp)
library(DRPT)
library(kerTests)

set.seed(231198)

## --- Shared parameters
n    = 150
m    = 150
MC   = 500
xxx  = 10
alphas = seq(0, 0.8, length.out = xxx)

gaussian.kernel = function(x, y, lambda = 1){
  d = length(x)
  lambda^(-d) * exp(-sum(((x - y)^2) / (lambda^2)))
}

rej.sampling.MMD = function(X, Y, S = 50, H = 99, r, M) {
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  n = nrow(X); dX = ncol(X); dY = ncol(Y)
  if (dX != dY) stop("X and Y must have same number of columns")
  
  if (!is.finite(M) || M <= 0) stop("M must be a positive finite envelope")
  
  ## vectorize r over rows of X
  r_vals = tryCatch({
    out = r(X[1, ])
    if (length(out) == 1L && is.finite(out)) apply(X, 1L, r) else as.numeric(r(X))
  }, error = function(e) apply(X, 1L, r))
  
  r_vals[!is.finite(r_vals)] = 0
  r_vals = pmax(r_vals, 0)
  if (any(r_vals > M + 1e-12)) {
    warning("Some r(x) > M; clipping to M to maintain validity.")
    r_vals = pmin(r_vals, M)
  }
  
  ## rejection step
  U = runif(n)
  keep = which(U <= (r_vals / M))
  if (length(keep) == 0L) return(1) 
  
  W = X[keep, , drop = FALSE]
  
  ## Two-sample test W vs Y via DRPT with r ≡ 1 (pure MMD + permutations)
  bw = kerTests::med_sigma(W, Y)
  kfun = function(u, v) gaussian.kernel(u, v, lambda = bw)
  
  r_one = function(x,y) as.numeric(0 * x + 1)
  
  DRPT::DRPT(W,Y, r = r_one, kernel = kfun)
}

## Study 1 (2d1): r(x,y) = 4 x y ; sup r = 4 on [0,1]^2
r1_xy  = function(x, y) 4 * x * y              
r1_vec = function(z)     4 * z[1] * z[2]      

invCDF = function(u) sqrt(u)

gen_Y_study1 = function(m, alpha){
  Y = matrix(0, nrow = m, ncol = 2)
  b = rbinom(m, 1, 1/(1 + alpha))
  n1 = sum(b == 1); n0 = m - n1
  if (n1 > 0) Y[b == 1, ] = cbind(invCDF(runif(n1)), invCDF(runif(n1)))
  if (n0 > 0) Y[b == 0, ] = cbind(rbeta(n0, 0.5, 0.5), rbeta(n0, 0.5, 0.5))
  Y
}

## Study 2 (2d2): r(x,y) = 2 x ; sup r = 2 on [0,1]^2
r2_xy  = function(x, y) 2 * x
r2_vec = function(z)     2 * z[1]

gen_Y_study2 = function(m, alpha){
  Y = matrix(0, nrow = m, ncol = 2)
  b = rbinom(m, 1, 1/(1 + alpha))
  n1 = sum(b == 1); n0 = m - n1
  if (n1 > 0) Y[b == 1, ] = cbind(invCDF(runif(n1)), runif(n1))
  if (n0 > 0) Y[b == 0, ] = cbind(rbeta(n0, 0.5, 0.5), rbeta(n0, 0.5, 0.5))
  Y
}

## ---------- Runner for one study ----------
run_study = function(study_id){
  pow_drpt = numeric(length(alphas))
  pow_rej  = numeric(length(alphas))
  decisions_drpt = vector("list", length(alphas))  
  decisions_rej  = vector("list", length(alphas))  
  
  for (i in seq_along(alphas)){
    alpha = alphas[i]
    cat(sprintf("[Study %d] alpha = %.3f\n", study_id, alpha))
    cnt_drpt = 0; cnt_rej = 0
    dec_drpt = integer(MC)  
    dec_rej  = integer(MC)
    
    for (b in 1:MC){
      ## Proposals X ~ f = Unif([0,1]^2)
      X = cbind(runif(n), runif(n))
      
      ## Generate Y
      if (study_id == 1){
        Y    = gen_Y_study1(m, alpha)
        r_xy = r1_xy
        r_v  = r1_vec
        Menv = 4
      } else {
        Y    = gen_Y_study2(m, alpha)
        r_xy = r2_xy
        r_v  = r2_vec
        Menv = 2
      }
      
      ## Kernel bandwidth (median heuristic on pooled)
      bw = kerTests::med_sigma(X, Y)
      kfun = function(u, v) gaussian.kernel(u, v, lambda = bw)
      
      ## DRPT
      p_drpt = DRPT::DRPT(X, Y, r = r_xy, kernel = kfun)
      
      ## Rejection-sampling MMD
      p_rej = rej.sampling.MMD(X = X, Y = Y, r = r_v, M = Menv)
      
      if (p_drpt < 0.05) { cnt_drpt = cnt_drpt + 1; dec_drpt[b] = 1 } else dec_drpt[b] = 0  
      if (p_rej  < 0.05) { cnt_rej  = cnt_rej  + 1; dec_rej[b]  = 1 } else dec_rej[b]  = 0  
    }
    
    pow_drpt[i] = cnt_drpt / MC
    pow_rej[i]  = cnt_rej  / MC
    decisions_drpt[[i]] = dec_drpt 
    decisions_rej[[i]]  = dec_rej 
  }
  
  list(alpha = alphas, power_drpt = pow_drpt, power_rej = pow_rej,
       decisions_drpt = decisions_drpt, decisions_rej = decisions_rej)
}

## ---------- Run both studies ----------
res1 = run_study(1)  # r=4xy
res2 = run_study(2)  # r=2x

## ---------- Save results ----------
dir.create("experiments/results", recursive = TRUE, showWarnings = FALSE)

write.csv(data.frame(Alpha = res1$alpha, Power = res1$power_drpt),
          "experiments/results/BIV1_DRPT.csv", row.names = FALSE)
write.csv(data.frame(Alpha = res1$alpha, Power_rej = res1$power_rej),
          "experiments/results/BIV1_REJ.csv",  row.names = FALSE)

write.csv(data.frame(Alpha = res2$alpha, Power = res2$power_drpt),
          "experiments/results/BIV2_DRPT.csv", row.names = FALSE)
write.csv(data.frame(Alpha = res2$alpha, Power_rej = res2$power_rej),
          "experiments/results/BIV2_REJ.csv",  row.names = FALSE)


## ---------- Save full 0/1 sequences (wide: one column per alpha, rows = MC iterations) ----------
alpha_cols = paste0("alpha_", format(alphas, nsmall = 3, trim = TRUE))

# Study 1
dec1_drpt = as.data.frame(do.call(cbind, res1$decisions_drpt)); names(dec1_drpt) = alpha_cols
dec1_rej  = as.data.frame(do.call(cbind, res1$decisions_rej));  names(dec1_rej)  = alpha_cols
dec1_drpt$iter = 1:MC; dec1_drpt = dec1_drpt[, c("iter", alpha_cols)]
dec1_rej$iter  = 1:MC; dec1_rej  = dec1_rej[,  c("iter", alpha_cols)]
write.csv(dec1_drpt, "experiments/results/BIV1_DRPT_decisions.csv", row.names = FALSE)
write.csv(dec1_rej,  "experiments/results/BIV1_REJ_decisions.csv",  row.names = FALSE)

# Study 2
dec2_drpt = as.data.frame(do.call(cbind, res2$decisions_drpt)); names(dec2_drpt) = alpha_cols
dec2_rej  = as.data.frame(do.call(cbind, res2$decisions_rej));  names(dec2_rej)  = alpha_cols
dec2_drpt$iter = 1:MC; dec2_drpt = dec2_drpt[, c("iter", alpha_cols)]
dec2_rej$iter  = 1:MC; dec2_rej  = dec2_rej[,  c("iter", alpha_cols)]
write.csv(dec2_drpt, "experiments/results/BIV2_DRPT_decisions.csv", row.names = FALSE)
write.csv(dec2_rej,  "experiments/results/BIV2_REJ_decisions.csv",  row.names = FALSE)
