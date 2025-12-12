rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")
# setwd("/storage/stats/strtng/SPT/synthetic")

library(MASS)
library(latex2exp)
library(DRPT)
library(kerTests)

set.seed(221198)

# Set parameters
n  = 150
m  = 150
d  = 1

MC = 500
S  = 50
xxx = 10
M  = 1

# Two shift functions
r4 = function(x) exp(-4 * x^2)
r1 = function(x) exp(-1 * x^2)

gaussian.kernel = function(x, y, lambda = 1){
  lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2)))
}

rej.sampling.MMD = function(X, Y, S = 50, H = 99, r, M) {
  # basic checks
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  n = nrow(X); m = nrow(Y)
  dX = ncol(X); dY = ncol(Y)
  if (dX != dY) stop("X and Y must have the same number of columns (dimension).")
  if (!is.finite(M) || M <= 0) stop("M must be a positive finite envelope with r(x) <= M.")
  
  # vectorize r over rows of X safely
  r_vals = tryCatch({
    out = r(X[1, ])               
    if (length(out) == 1L && is.finite(out)) {
      # apply row-wise
      apply(X, 1L, r)
    } else {
      # user-supplied r already vectorizes? try calling on full matrix
      as.numeric(r(X))
    }
  }, error = function(e) {
    # fallback to row-wise apply
    apply(X, 1L, r)
  })
  
  # rejection sampling (thinning)
  U = runif(n)
  keep = which(U <= (r_vals / M))
  if (length(keep) == 0L) {
    return(1)
  }
  W = X[keep, , drop = FALSE]
  
  ## Two-sample test W vs Y via DRPT with r ≡ 1 (pure MMD + permutations)
  bw = kerTests::med_sigma(W, Y)
  kfun = function(u, v) gaussian.kernel(u, v, lambda = bw)
  
  r_one = function(x) as.numeric(0 * x + 1)
  
  # Two-sample test W vs Y (null: same distribution) using your DRPT with r ≡ 1
  p_hat = DRPT::DRPT(W, Y, r_one, kernel = kfun)
  
  return(p_hat)
}


# Prepare containers
alphas = seq(from = 0, to = 0.4, length.out = xxx)

pow_drpt_r4    = numeric(xxx)
pow_drpt_r1    = numeric(xxx)
pow_rej_r4     = numeric(xxx)
pow_rej_r1     = numeric(xxx)

alpha_index = 1

# Loop over alpha values
for (alpha in alphas){
  message(sprintf("alpha = %.4f", alpha))
  
  # Counters for rejections
  cnt_drpt_r4 = 0
  cnt_drpt_r1 = 0
  cnt_rej_r4  = 0
  cnt_rej_r1  = 0
  
  # Monte Carlo
  for (b in 1:MC) {
    # X ~ N(0,1)
    X =as.matrix(rnorm(n, 0, 1))
    
    # Y: mixture; vectorized draw of component labels for r4
    comp = rbinom(m, 1, 1/(1 + alpha))  # 1 => Gaussian component
    Yvec = ifelse(comp == 1, rnorm(m, 0, 1/3), rexp(m, 1))
    Y4 = as.matrix(Yvec)
    
    # Y: mixture; vectorized draw of component labels for r1
    comp = rbinom(m, 1, 1/(1 + alpha))  # 1 => Gaussian component
    Yvec = ifelse(comp == 1, rnorm(m, 0, sqrt(1/3)), rexp(m, 1))
    Y1 = as.matrix(Yvec)
    
    # compute median bandwidth
    median.bandwidth4 = kerTests::med_sigma(X, Y4)
    kfun4 = function(x, y) gaussian.kernel(x, y, lambda = median.bandwidth4)
    
    median.bandwidth1 = kerTests::med_sigma(X, Y1)
    kfun1 = function(x, y) gaussian.kernel(x, y, lambda = median.bandwidth1)
    
    # --- run the tests ---
    p_drpt_r4 = DRPT(X, Y4, r = r4, kernel = kfun4)
    p_drpt_r1 = DRPT(X, Y1, r = r1, kernel = kfun1)
    
    p_rej_r4 = rej.sampling.MMD(X, Y4, r = r4, M = M)
    p_rej_r1 = rej.sampling.MMD(X, Y1, r = r1, M = M)
    
    if (p_drpt_r4 < 0.05) cnt_drpt_r4 =cnt_drpt_r4 + 1
    if (p_drpt_r1 < 0.05) cnt_drpt_r1 =cnt_drpt_r1 + 1
    if (p_rej_r4 < 0.05) cnt_rej_r4 =cnt_rej_r4 + 1
    if (p_rej_r1 < 0.05) cnt_rej_r1 =cnt_rej_r1 + 1
  }
  
  # Power estimates at this alpha
  pow_drpt_r4[alpha_index] = cnt_drpt_r4 / MC
  pow_drpt_r1[alpha_index] = cnt_drpt_r1 / MC
  pow_rej_r4[alpha_index]  = cnt_rej_r4  / MC
  pow_rej_r1[alpha_index]  = cnt_rej_r1  / MC
  
  print(cnt_drpt_r4 / MC)
  print(cnt_drpt_r1 / MC)
  print(cnt_rej_r4  / MC)
  print(cnt_rej_r1  / MC)
  
  alpha_index =alpha_index + 1
}

# --- Save results to the four files your plotting code expects ---

# DRPT - exp(-4x^2)
dir.create("experiments/results", recursive = TRUE, showWarnings = FALSE)
write.csv(
  data.frame(Alpha = alphas, Power = pow_drpt_r4),
  "experiments/results/simul_1D.csv",
  row.names = FALSE
)

# DRPT - exp(-x^2)
write.csv(
  data.frame(Alpha = alphas, Power = pow_drpt_r1),
  "experiments/results/simul1D_x2.csv",
  row.names = FALSE
)

# MMD Rejection - exp(-4x^2)
write.csv(
  data.frame(Alpha = alphas, Power_rej = pow_rej_r4),
  "experiments/results/simul_rej1D.csv",
  row.names = FALSE
)

# MMD Rejection - exp(-x^2)
write.csv(
  data.frame(Alpha = alphas, Power_rej = pow_rej_r1),
  "experiments/results/simul_rej1D_e-x^2.csv",
  row.names = FALSE
)

# --- Plot (unchanged inputs) ---
data.MMD  = read.csv(file = "experiments/results/simul_1D.csv")
data.x2   = read.csv(file = "experiments/results/simul1D_x2.csv")
data.rej  = read.csv(file = "experiments/results/simul_rej1D.csv")
data.rej2 = read.csv(file = "experiments/results/simul_rej1D_e-x^2.csv")

# MC must be defined in your environment (number of Monte Carlo runs)
# e.g., MC = 1000

compute_sd = function(power, alpha_vec) {
  m = power
  sqrt(m * (1 - m) / MC)
}

sd1_drpt = compute_sd(data.MMD$Power,  data.MMD$Alpha)
sd1_rej  = compute_sd(data.rej$Power_rej,  data.MMD$Alpha)
sd2_drpt = compute_sd(data.x2$Power,   data.MMD$Alpha)
sd2_rej  = compute_sd(data.rej2$Power_rej, data.MMD$Alpha)

# Helper to draw error bars, clipped to [0,1]
add_err = function(x, y, sd, col) {
  ylow  = pmax(0, y - sd)
  yhigh = pmin(1, y + sd)
  arrows(x0 = x, y0 = ylow, x1 = x, y1 = yhigh,
         angle = 90, code = 3, length = 0.04,
         col = col, lwd = 1.3)
}

# Compute y-limits considering error bars
all_low  = pmin(data.MMD$Power  - sd1_drpt,
                 data.x2$Power   - sd2_drpt,
                 data.rej$Power_rej  - sd1_rej,
                 data.rej2$Power_rej - sd2_rej, na.rm = TRUE)
all_high = pmax(data.MMD$Power  + sd1_drpt,
                 data.x2$Power   + sd2_drpt,
                 data.rej$Power_rej  + sd1_rej,
                 data.rej2$Power_rej + sd2_rej, na.rm = TRUE)
ylim_use = c(max(0, min(all_low, na.rm = TRUE)),
              min(1, max(all_high, na.rm = TRUE)))

dir.create("experiments/pictures", recursive = TRUE, showWarnings = FALSE)
png("experiments/pictures/simul_1D.png")
par(mfrow = c(1, 1))

plot(data.MMD$Alpha, data.MMD$Power, type = "b",
     col = "darkviolet", lwd = 2, pch = 14,
     xlab = TeX(r'($\eta$)', bold = TRUE), ylab = "Power",
     ylim = ylim_use)

# Other series
lines(data.x2$Alpha,  data.x2$Power,          type = "b", col = "orange", lwd = 2, pch = 14)
lines(data.rej$Alpha, data.rej$Power_rej,     type = "b", col = "green",  lwd = 2, pch = 16)
lines(data.rej2$Alpha,data.rej2$Power_rej,    type = "b", col = "blue",   lwd = 2, pch = 16)

# SD segment bars
add_err(data.MMD$Alpha,  data.MMD$Power,        sd1_drpt, "darkviolet")
add_err(data.x2$Alpha,   data.x2$Power,         sd2_drpt, "orange")
add_err(data.rej$Alpha,  data.rej$Power_rej,    sd1_rej,  "green")
add_err(data.rej2$Alpha, data.rej2$Power_rej,   sd2_rej,  "blue")

legend("topleft", inset = 0.01, horiz = FALSE, lty = 1, bty = "n",
       legend = c("(E1) - exp(-4x^2)", "(E1) - exp(-x^2)",
                  "(E2) - exp(-4x^2)", "(E2) - exp(-x^2)"),
       col = c("darkviolet", "orange", "green", "blue"),
       pch = c(14, 14, 16, 16))

abline(h = 0.05, col = "red", lty = 2)
dev.off()