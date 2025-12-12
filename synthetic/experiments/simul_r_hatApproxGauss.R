rm(list = ls())  # Clear environment
gc()             # Free memory

# setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")
setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/synthetic/")

library(MASS)
library(latex2exp)
library(squash)
library(DRPT)
library(kerTests)

set.seed(110932)

# Set parameters
n = 250
m = 250
d = 1

MC = 300
S = 50

gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

invCDF = function(x){
  return(sqrt(x))
}


#---------------------------------
# Estimate LL
#---------------------------------
size_n = numeric(4)
n_index = 1
alpha = 0

decisions_all = list()

for (eps in c(0.05, 0.10, 0.15, 0.20, 0.25)){
  
  r_N  = function(x) return(exp((1+eps)*x))
  
  sum_ind = 0
  decisions_seq = integer(MC)
  
  #--------------------------------------
  # Monte Carlo simulation on test data
  #--------------------------------------
  for (b in 1:MC) {
    ######## generate data
    X = rnorm(n, 0, 1)
    Y = rnorm(n, 1, 1)
    
    bw = kerTests::med_sigma(X, Y)
    kfun = function(u, v) gaussian.kernel(u, v, lambda = bw)
    
    ########### Our test: call the SPT function
    p_M_alpha = DRPT(X, Y, r = r_N, kernel = kfun)
    
    if (p_M_alpha < 0.05) {
      sum_ind = sum_ind + 1
      decisions_seq[b] = 1      
    } else {
      decisions_seq[b] = 0     
    }
  }
  
  # Save the whole 0/1 sequence for this eps
  decisions_all[[paste0("eps_", sprintf("%.2f", eps))]] = decisions_seq
  
  # Calculate the power for this alpha
  size_n[n_index] = sum_ind / MC
  print(sum_ind / MC)
  n_index = n_index + 1
}



################# PLOT #########################
# Load saved results (means and 0/1 sequences)
means_df = read.csv("experiments/results/r_estimation_r_N_Gaussian.csv")

# ---- Parameters & helpers ----
delta  = c(0.05, 0.10, 0.15, 0.20, 0.25)
nu     = 1 + delta
theory = 0.05 + c(0.31, 0.57, 0.76, 0.89, 0.95)

# Mean power from simulations (match deltas to eps in means_df$n)
idx  = match(delta, means_df$n)
drpt = means_df$Size[idx]

# SD from means: sqrt(mean * (1 - mean) / MC)
sd_drpt = sqrt(drpt * (1 - drpt) / MC)

# Order so the largest Delta is at the top
ord = order(delta, decreasing = TRUE)

# Ensure output dir exists (MATCHES path below)
dir.create("experiments/pictures", recursive = TRUE, showWarnings = FALSE)

png(filename = "experiments/pictures/DRPT_excessGaussian.png")
op = par(mar = c(4, 9, 3, 1))

y     = seq_along(delta)
ylabs = paste0("\u03bc + ", sprintf("%.2f", delta[ord]))

# empty canvas
plot(NA, xlim = c(0, 1), ylim = c(0.5, length(y) + 0.5),
     xlab = "Type-I error", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

# axes & grid
axis(1, at = seq(0, 1, 0.1), labels = sprintf("%.1f", seq(0, 1, 0.1)))
axis(2, at = y, labels = ylabs, las = 1, tick = FALSE)
abline(h = y, col = "grey90", lwd = 1)

# connectors
for (i in seq_along(y)) {
  segments(min(drpt[ord][i], theory[ord][i]), y[i],
           max(drpt[ord][i], theory[ord][i]), y[i],
           col = "grey70", lwd = 3)
}

# points
points(theory[ord], y, pch = 16, cex = 1.2, col = "red")
points(drpt[ord],   y, pch = 16, cex = 1.2, col = "purple")

# --- SD lines for DRPT (horizontal at each y), clipped to [0,1] ---
for (i in seq_along(y)) {
  x0 = max(0, drpt[ord][i] - sd_drpt[ord][i])
  x1 = min(1, drpt[ord][i] + sd_drpt[ord][i])
  segments(x0, y[i], x1, y[i], col = "purple", lwd = 3)
}

# reference line at alpha=0.05
abline(v = 0.05, lty = 2)

# labels & legend
mtext(expression(nu == mu + Delta), side = 2, line = 6.5, las = 0)
legend("topright", inset = 0.02, bty = "n",
       legend = c("Theoretical bound", "DRPT (E1)"),
       pch = 16, col = c("red", "purple"), pt.cex = 1.2)

par(op)
dev.off()

# # --------------------------------------
# # combine pictures
# # --------------------------------------
# 
# library(magick)
# 
# files = c(
#   "experiments/pictures/BIVsimul_shiftedMMD.png",
#   "experiments/pictures/binary.png",
#   "experiments/pictures/DRPT_excessGaussian.png"
# )
# imgs = image_read(files)
# 
# # Draw (a)/(b)/(c)
# label_draw = function(img, lab, x = 15, y = 35) {
#   info = image_info(img)
#   h = info$height
#   cex_val = max(1, h / 400) 
#   
#   img = image_draw(img)  
#   par(mar = c(0,0,0,0), xpd = NA, family = "HersheySans")
#   # soft white box behind text
#   w  = strwidth(lab, cex = cex_val)
#   ht = strheight(lab, cex = cex_val)
#   rect(x - 6, y - 1.6*ht, x + 6 + 1.1*w, y + 6,
#        col = rgb(1,1,1,0.7), border = NA)
#   text(x, y, labels = lab, cex = cex_val, col = "black", adj = c(0, 1))
#   dev.off()
#   img  
# }
# 
# a = label_draw(imgs[1], "(a)")
# b = label_draw(imgs[2], "(b)")
# c = label_draw(imgs[3], "(c)")
# 
# # Combine horizontally
# combo = image_append(image_join(a, b, c), stack = FALSE)
# image_write(combo, "experiments/pictures/combined_1.png")
