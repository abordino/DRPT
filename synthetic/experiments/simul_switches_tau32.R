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

N  = 100
MC = 2000

eta_grid = c(0, 0.9)
tau_grid = c(0.1, 0.3, 0.5, 0.7, 0.9)

H = 99

# keep the number of proposed pairs per Markov step fixed at 5,
# and vary the number of Markov steps S from 1 to 10.
number_pairs_switched = 5
S_grid = 1:10

# Total number of pair-switch proposals across the S Markov steps.
# This is the quantity shown on the x-axis.
number_switches_grid = number_pairs_switched * S_grid

level = 0.05


## ---------- Kernel ----------

gaussian.kernel = function(x, y, lambda = 1){
  d = length(x)
  lambda^(-d) * exp(-sum(((x - y)^2) / (lambda^2)))
}


## ---------- Study 1: r(x,y) = 4xy ----------

r_xy = function(x, y){
  4 * x * y
}

invCDF = function(u){
  sqrt(u)
}


## ---------- Y generator ----------
## eta = 0 is the null; eta = 0.9 is the alternative.
## The alternative component is Beta(1.5, 1.5) in each coordinate.

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
    length(S_grid),
    length(tau_grid),
    length(eta_grid)
  ),
  dimnames = list(
    paste0("S_", S_grid),
    paste0("tau_", tau_grid),
    paste0("eta_", eta_grid)
  )
)


decisions = array(
  0,
  dim = c(
    MC,
    length(S_grid),
    length(tau_grid),
    length(eta_grid)
  ),
  dimnames = list(
    paste0("iter_", 1:MC),
    paste0("S_", S_grid),
    paste0("tau_", tau_grid),
    paste0("eta_", eta_grid)
  )
)


## ---------- Run simulation ----------

for (a in seq_along(eta_grid)){
  
  eta = eta_grid[a]
  
  for (t in seq_along(tau_grid)){
    
    tau = tau_grid[t]
    
    n = round(N * tau)
    m = N - n
    
    for (s_idx in seq_along(S_grid)){
      
      S = S_grid[s_idx]
      total_pair_proposals = number_pairs_switched * S
      
      cat(
        sprintf(
          paste0(
            "eta = %.1f | tau = %.1f | n = %d | m = %d | ",
            "H = %d | S = %d | pairs/step = %d | total proposals = %d\n"
          ),
          eta, tau, n, m, H, S,
          number_pairs_switched,
          total_pair_proposals
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
          S = S,
          number_pairs_switched = number_pairs_switched
        )
        
        dec[b] = as.integer(
          p_drpt <= level
        )
      }
      
      decisions[, s_idx, t, a] = dec
      power[s_idx, t, a] = mean(dec)
      
      cat(
        sprintf(
          "Power = %.3f\n\n",
          power[s_idx, t, a]
        )
      )
    }
  }
}


## ---------- Results ----------

results_list = list()
rr = 1

for (a in seq_along(eta_grid)){
  
  for (t in seq_along(tau_grid)){
    
    tau = tau_grid[t]
    n = round(N * tau)
    m = N - n
    
    results_list[[rr]] = data.frame(
      Eta = eta_grid[a],
      Tau = tau,
      n = n,
      m = m,
      H = H,
      S = S_grid,
      Number_pairs_per_step = number_pairs_switched,
      Total_pair_proposals = number_switches_grid,
      Power = power[, t, a]
    )
    
    rr = rr + 1
  }
}

results = do.call(
  rbind,
  results_list
)

rownames(results) = NULL


## ---------- Monte Carlo standard errors ----------

results$SE = sqrt(
  results$Power *
    (1 - results$Power) /
    MC
)

print(results)


## ---------- Save results ----------

dir.create(
  "experiments/results",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  results,
  "experiments/results/BIV1_DRPT_steps5_tau_eta09_32.csv",
  row.names = FALSE
)


## ---------- Save full 0/1 decisions ----------

decisions_df = data.frame(
  iter = 1:MC
)

for (a in seq_along(eta_grid)){
  
  for (t in seq_along(tau_grid)){
    
    tau = tau_grid[t]
    n = round(N * tau)
    m = N - n
    
    for (s_idx in seq_along(S_grid)){
      
      S = S_grid[s_idx]
      total_pair_proposals = number_pairs_switched * S
      
      col_name = paste0(
        "eta_", eta_grid[a],
        "_tau_", tau,
        "_n_", n,
        "_m_", m,
        "_H_", H,
        "_S_", S,
        "_pairsPerStep_", number_pairs_switched,
        "_totalProposals_", total_pair_proposals
      )
      
      decisions_df[[col_name]] =
        decisions[, s_idx, t, a]
    }
  }
}

write.csv(
  decisions_df,
  "experiments/results/BIV1_DRPT_steps5_tau_eta09_32_decisions.csv",
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
    lwd = 1.2
  )
}


## ---------- Plot power vs total number of pair-switch proposals ----------

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
  "experiments/pictures/BIV1_DRPT_power_vs_steps5_tau_eta09_32.png",
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


## ---------- Empty plotting region ----------

plot(
  NA,
  xlim = range(number_switches_grid),
  ylim = c(0, 1),
  xaxt = "n",
  xlab = "Number of proposed pairwise switches",
  ylab = "Power",
  cex.lab = 1.25,
  cex.axis = 1.1,
  cex.main = 1.3
)

axis(
  1,
  at = number_switches_grid,
  labels = number_switches_grid,
  cex.axis = 1.05
)


## ---------- Add all curves ----------

for (a in seq_along(eta_grid)){
  
  for (t in seq_along(tau_grid)){
    
    xvals = number_switches_grid
    yvals = power[, t, a]
    sevals = power_se[, t, a]
    
    lines(
      xvals,
      yvals,
      type = "b",
      col = cols[t],
      lwd = 2.5,
      lty = ltys[a],
      pch = pchs[t],
      cex = 1.15
    )
    
    add_err(
      x = xvals,
      y = yvals,
      se = sevals,
      col = cols[t]
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

## ---------- Highlight 5*S = min(n,m), eta = 0.9 ----------

a = which(eta_grid == 0.9)

for (t in seq_along(tau_grid)) {
  
  n = round(N * tau_grid[t])
  m = N - n
  
  x_target = min(n, m)
  k = which(5 * S_grid == x_target)
  
  points(
    x_target,
    power[k, t, a],
    pch = 1,
    col = "red",
    cex = 2,
    lwd = 2
  )
}


## ---------- Legends ----------

legend(
  "topleft",
  inset = 0.01,
  legend = parse(
    text = paste0("n/(n+m) == ", tau_grid)
  ),
  col = cols,
  pch = pchs,
  lty = 1,
  lwd = 2.5,
  pt.cex = 1.15,
  cex = 0.95,
  bty = "n",
  title = expression(n/(n + m))
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
