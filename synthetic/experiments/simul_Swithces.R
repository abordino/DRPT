rm(list = ls())
gc()

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")
# setwd("/storage/stats/strtng/SPT/synthetic")

library(DRPT)
library(kerTests)

source("experiments/DRPT_v2.R")
source("experiments/DRPT_v3.R")

set.seed(231198)


## ---------- Parameters ----------

n  = 50
m  = 50
MC = 2000

eta_grid = c(0, 0.9)

level = 0.05


S_fixed = 1
number_pairs_grid = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
H_grid = c(39, 59, 79, 99, 119)


## ---------- Kernel ----------

gaussian.kernel = function(x, y, lambda = 1){
  
  d = length(x)
  
  lambda^(-d) *
    exp(
      -sum(
        ((x - y)^2) / (lambda^2)
      )
    )
}


## ---------- Study 1: r(x,y) = 4xy ----------

r_xy = function(x, y){
  4 * x * y
}


invCDF = function(u){
  sqrt(u)
}


gen_Y = function(m, eta){
  
  Y = matrix(
    0,
    nrow = m,
    ncol = 2
  )
  
  b = rbinom(
    m,
    1,
    1 / (1 + eta)
  )
  
  n1 = sum(b == 1)
  n0 = m - n1
  
  
  if (n1 > 0){
    
    Y[b == 1, ] = cbind(
      invCDF(runif(n1)),
      invCDF(runif(n1))
    )
  }
  
  
  if (n0 > 0){
    
    Y[b == 0, ] = cbind(
      rbeta(n0, 1.5, 1.5),
      rbeta(n0, 1.5, 1.5)
    )
  }
  
  
  Y
}


## ---------- Storage ----------

power = array(
  0,
  dim = c(
    length(number_pairs_grid),
    length(H_grid),
    length(eta_grid)
  ),
  dimnames = list(
    paste0("pairs_", number_pairs_grid),
    paste0("H_", H_grid),
    paste0("eta_", eta_grid)
  )
)


decisions = array(
  0,
  dim = c(
    MC,
    length(number_pairs_grid),
    length(H_grid),
    length(eta_grid)
  ),
  dimnames = list(
    paste0("iter_", 1:MC),
    paste0("pairs_", number_pairs_grid),
    paste0("H_", H_grid),
    paste0("eta_", eta_grid)
  )
)


## ---------- Run simulation ----------

for (a in seq_along(eta_grid)){
  
  eta = eta_grid[a]
  
  
  for (h in seq_along(H_grid)){
    
    H = H_grid[h]
    
    
    for (k in seq_along(number_pairs_grid)){
      
      number_pairs_switched = number_pairs_grid[k]
      
      
      cat(
        sprintf(
          paste0(
            "eta = %.3f | H = %d | S = %d | ",
            "number_pairs_switched = %d\n"
          ),
          eta,
          H,
          S_fixed,
          number_pairs_switched
        )
      )
      
      
      dec = integer(MC)
      
      
      for (b in 1:MC){
        
        X = cbind(
          runif(n),
          runif(n)
        )
        
        
        Y = gen_Y(
          m,
          eta
        )
        
        
        bw = kerTests::med_sigma(
          X,
          Y
        )
        
        
        kfun = function(u, v){
          
          gaussian.kernel(
            u,
            v,
            lambda = bw
          )
        }
        
        
        p_drpt = DRPT_v3(
          X,
          Y,
          r = r_xy,
          kernel = kfun,
          H = H,
          S = S_fixed,
          number_pairs_switched = number_pairs_switched
        )
        
        
        dec[b] = as.integer(
          p_drpt <= level
        )
      }
      
      
      decisions[, k, h, a] = dec
      
      power[k, h, a] = mean(dec)
      
      
      cat(
        sprintf(
          "Power = %.3f\n\n",
          power[k, h, a]
        )
      )
    }
  }
}


## ---------- Results ----------

results = do.call(
  rbind,
  lapply(
    seq_along(eta_grid),
    function(a){
      
      do.call(
        rbind,
        lapply(
          seq_along(H_grid),
          function(h){
            
            data.frame(
              Eta = eta_grid[a],
              H = H_grid[h],
              S = S_fixed,
              Number_pairs_switched = number_pairs_grid,
              Power = power[, h, a]
            )
          }
        )
      )
    }
  )
)


rownames(results) = NULL


## ---------- Monte Carlo standard errors ----------

results$SE = sqrt(
  results$Power *
    (1 - results$Power) /
    MC
)


## ---------- Print results ----------

print(results)


## ---------- Save results ----------

dir.create(
  "experiments/results",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  results,
  "experiments/results/BIV1_DRPT_pairs_H_eta_S1.csv",
  row.names = FALSE
)


## ---------- Save full 0/1 decisions ----------

decisions_df = data.frame(
  iter = 1:MC
)


for (a in seq_along(eta_grid)){
  
  for (h in seq_along(H_grid)){
    
    for (k in seq_along(number_pairs_grid)){
      
      col_name = paste0(
        "eta_", eta_grid[a],
        "_H_", H_grid[h],
        "_S_", S_fixed,
        "_pairs_", number_pairs_grid[k]
      )
      
      
      decisions_df[[col_name]] =
        decisions[, k, h, a]
    }
  }
}


write.csv(
  decisions_df,
  "experiments/results/BIV1_DRPT_pairs_H_eta_S1_decisions.csv",
  row.names = FALSE
)


## ---------- Standard errors ----------

power_se = sqrt(
  power * (1 - power) / MC
)

add_err = function(x, y, se, col){
  
  ylow = pmax(
    0,
    y - se
  )
  
  yhigh = pmin(
    1,
    y + se
  )
  
  arrows(
    x0 = x,
    y0 = ylow,
    x1 = x,
    y1 = yhigh,
    angle = 90,
    code = 3,
    length = 0.04,
    col = col,
    lwd = 1.4
  )
}


## ---------- Plot power as a function of number of pairs switched ----------

dir.create(
  "experiments/pictures",
  recursive = TRUE,
  showWarnings = FALSE
)


cols = c(
  "darkviolet",
  "green",
  "blue",
  "brown", 
  "orange"
)

pchs = c(
  14,
  16,
  17,
  18, 
  19
)


ltys = c(
  2,
  1
)


## ---------- PNG ----------

png(
  "experiments/pictures/BIV1_DRPT_power_vs_pairs_eta_S1.png",
  width = 8,
  height = 6.5,
  units = "in",
  res = 400,
  pointsize = 12
)


par(
  mfrow = c(1, 1),
  mar = c(5, 5, 4, 2) + 0.1,
  las = 1
)


## ---------- Initial curve ----------

plot(
  number_pairs_grid,
  power[, 1, 1],
  type = "b",
  col = cols[1],
  lwd = 2.5,
  lty = ltys[1],
  pch = pchs[1],
  cex = 1.15,
  ylim = c(0, 1),
  xaxt = "n",
  xlab = "Number of proposed pairwise switches",
  ylab = "Power",
  cex.lab = 1.25,
  cex.axis = 1.1,
  cex.main = 1.3
)


## ---------- X axis ----------

axis(
  1,
  at = number_pairs_grid,
  labels = number_pairs_grid,
  cex.axis = 1.05
)


## ---------- Add all remaining curves ----------

for (a in seq_along(eta_grid)){
  
  for (h in seq_along(H_grid)){
    
    if (a == 1 && h == 1){
      next
    }
    
    
    lines(
      number_pairs_grid,
      power[, h, a],
      type = "b",
      col = cols[h],
      lwd = 2.5,
      lty = ltys[a],
      pch = pchs[h],
      cex = 1.15
    )
  }
}


## ---------- Add Monte Carlo standard-error bars ----------

for (a in seq_along(eta_grid)){
  
  for (h in seq_along(H_grid)){
    
    add_err(
      x = number_pairs_grid,
      y = power[, h, a],
      se = power_se[, h, a],
      col = cols[h]
    )
  }
}


## ---------- Significance level ----------

abline(
  h = level,
  col = "red",
  lty = 2,
  lwd = 1.5
)


## ---------- Legends ----------

legend(
  "topleft",
  inset = 0.01,
  legend = paste0("H = ", H_grid),
  col = cols,
  pch = pchs,
  lty = 1,
  lwd = 2.5,
  pt.cex = 1.15,
  cex = 1.0,
  bty = "n",
  title = "H"
)


legend(
  "topright",
  inset = 0.01,
  legend = c(
    expression(eta == 0),
    expression(eta == 0.9)
  ),
  col = "black",
  lty = c(2, 1),
  lwd = 2.5,
  pch = NA,
  cex = 1.0,
  bty = "n",
  title = expression(eta)
)


dev.off()