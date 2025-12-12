rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/twoMoon")
# setwd("/storage/stats/strtng/SPT/twoMoon")

library(MASS)
library(latex2exp)
library(DRPT)
library(kerTests)
library(reticulate)

#--------------------------------------------------------------------
# Want to run it three times with these seeds, then combine
set.seed(170966)
# set.seed(180966)
# set.seed(200801)
#--------------------------------------------------------------------

# DRPT test using a TorchScript r_hat(theta, x) exported by train_moon_lampe.py

# 1) Load TorchScript model via reticulate/torch -------------------------------
virtualenv_create("sbi-env")
use_virtualenv("sbi-env", required = TRUE)
system('pip install zuko lampe torch')
pt = import("torch")

rhat_loader = function(path) {
  mod = pt$jit$load(path)
  mod$eval()
  
  # helper: broadcast two numeric vectors/scalars to length n
  bcast2 = function(a, b, n) cbind(rep(a, length.out = n), rep(b, length.out = n))
  
  function(theta1, theta2, x1, x2) {
    # coerce to numeric
    theta1 = as.numeric(theta1); theta2 = as.numeric(theta2)
    x1     = as.numeric(x1);     x2     = as.numeric(x2)
    
    # determine common length
    lens = c(length(theta1), length(theta2), length(x1), length(x2))
    n = max(lens)
    if (!all(lens %in% c(1, n))) {
      stop("All inputs must be scalars or length-n vectors with the same n.")
    }
    
    # build n×2 matrices
    theta = bcast2(theta1, theta2, n)
    x     = bcast2(x1, x2, n)
    
    th_t = pt$tensor(theta, dtype = pt$float32)
    xx_t = pt$tensor(x,     dtype = pt$float32)
    
    # no-grad + CPU + NumPy
    old = pt$is_grad_enabled()
    pt$set_grad_enabled(FALSE); on.exit(pt$set_grad_enabled(old), add = TRUE)
    out = mod$forward(th_t, xx_t)
    as.numeric(out$detach()$cpu()$numpy())
  }
}

rsim_theta = function(n) cbind(runif(n, -1, 1), runif(n, -1, 1))
rsim_x_given_theta = function(theta) {
  n = nrow(theta)
  alpha = runif(n, -pi/2, pi/2)
  rr = rnorm(n, mean = 0.1, sd = 0.01)
  base  = cbind(rr * cos(alpha) + 0.25, rr * sin(alpha))
  shift = cbind(-abs(theta[,1] + theta[,2]) / sqrt(2), (-theta[,1] + theta[,2]) / sqrt(2))
  base + shift
}

# Choose estimator: BNRE or NRE
N = 200
d = 4
gaussian.kernel = function(x, y, lambda = 1){
  lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2)))
}
# ---- main evaluation over budgets & methods ----
budgets = c(2^4, 2^7, 2^9, 2^12, 2^15)
methods = c("bnre", "nre")  

reps = 100

dir.create("results", showWarnings = FALSE, recursive = TRUE)

all_res = list()

start = Sys.time()
for (B in budgets) {
  print(B)
  for (m in methods) {
    print(m)
    
    model_path = sprintf("artifacts/moon_rhat_%s_B%d.pt", m, B)
    if (!file.exists(model_path)) {
      warning(sprintf("Missing model file: %s — skipping.", model_path))
      next
    }
    
    r_hat = rhat_loader(model_path)
    
    create_lookup_list = function(X, rX, Y, rY) {
      keysX = apply(X, 1, paste, collapse = "_")
      keysY = apply(Y, 1, paste, collapse = "_")
      
      lookup_list = setNames(as.list(c(rX, rY)), c(keysX, keysY))
      
      return(lookup_list)
    }
    
    pvals = numeric(reps)
    rejs  = logical(reps)
    
    for (b in 1:reps){
      
      thJ = rsim_theta(N); xJ = rsim_x_given_theta(thJ)              # joint
      thP = rsim_theta(N); xP = rsim_x_given_theta(rsim_theta(N))    # product-of-marginals
      
      Y = cbind(thJ, xJ)   # (theta, x) from joint
      X = cbind(thP, xP)   # (theta, x) from product
      
      Z = rbind(X,Y)
      r.values = r_hat(Z[,1], Z[,2], Z[,3], Z[,4])
      lookup_table = create_lookup_list(X, r.values[1:N], Y, tail(r.values, N))
      
      r.hat = function(...) {
        z = cbind(...)
        if (is.matrix(z) || is.data.frame(z)) {
          keys = apply(z, 1, paste, collapse = "_")
          return(sapply(keys, function(key) ifelse(key %in% names(lookup_table), lookup_table[[key]], NA)))
        } else {
          key = paste(z, collapse = "_")
          return(ifelse(key %in% names(lookup_table), lookup_table[[key]], NA))
        }
      }
      
      median.bandwidth = kerTests::med_sigma(X, Y)
      kfun = function(x, y) gaussian.kernel(x, y, lambda = median.bandwidth)
      
      # run DRPT
      res = DRPT(X = X, Y = Y, r = r.hat, kernel = kfun, H = 99, S = 50)
      
      # p-value extraction (numeric or list with p.value)
      p = if (is.list(res) && !is.null(res$p.value)) res$p.value else as.numeric(res)
      if (!is.finite(p)) p = NA_real_
      
      pvals[b] = p
      rejs[b]  = is.finite(p) && (p < 0.05)
    }
    
    # Save raw replicate-level results (0/1 rejections + p-values)
    raw_df = data.frame(
      method = toupper(m),
      budget = B,
      rep    = seq_len(reps),
      pval   = pvals,
      reject = as.integer(rejs)  # 0/1 sequence
    )
    raw_outfile <- "results/drpt_trials_1.csv"
    utils::write.table(
      raw_df,
      file       = raw_outfile,
      sep        = ",",
      row.names  = FALSE,
      col.names  = !file.exists(raw_outfile),
      append     = file.exists(raw_outfile)
    )
  }
}

#--------------------------------------------------------------------
# Run the simulations with the three different seeds
# Then load the results and combine them
#--------------------------------------------------------------------

budgets = c(2^4, 2^7, 2^9, 2^12, 2^15)
methods = c("bnre", "nre")  

reps = 300 # now have triple the sample size

suppressPackageStartupMessages({ library(ggplot2) })

all_trials_1 = utils::read.csv("results/drpt_trials_1.csv", header = FALSE)
all_trials_2 = utils::read.csv("results/drpt_trials_2.csv", header = FALSE)
all_trials_3 = utils::read.csv("results/drpt_trials_3.csv", header = FALSE)

all_trials = rbind(all_trials_1, all_trials_2, all_trials_3)

all_res = list()
for (B in budgets) {
  for (m in methods) {
    rejs = all_trials[tolower(all_trials$V1) == m & all_trials$V2 == B, ]$V5
    # Per-(method,budget) summaries with SDs
    all_res[[length(all_res) + 1]] = data.frame(
      method = toupper(m),
      budget = B,
      power  = mean(rejs,  na.rm = TRUE),                 # mean of 0/1 = power
      sd_rej = sqrt(mean(rejs,  na.rm = TRUE)*(1-mean(rejs,  na.rm = TRUE))/reps)         # SD of 0/1 sequence
    )
  }
}

summary_df = do.call(rbind, all_res)
utils::write.csv(summary_df, file = "results/drpt_summary.csv", row.names = FALSE)


summary_df = utils::read.csv("results/drpt_summary.csv")

# clamp error bars to [0,1] ranges where appropriate
summary_df$power_lo = pmax(0, summary_df$power - summary_df$sd_rej)
summary_df$power_hi = pmin(1, summary_df$power + summary_df$sd_rej)

# dodge so methods at the same budget don't sit exactly on top of each other
pos = ggplot2::position_dodge(width = 0.08)

p_power = ggplot2::ggplot(
  summary_df,
  ggplot2::aes(x = budget, y = power, color = method, group = method)
) +
  ggplot2::geom_line(position = pos, linewidth = 1) +
  ggplot2::geom_point(
    ggplot2::aes(shape = method),
    size = 2.7,
    position = pos,
    stroke = 0.6
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = power_lo, ymax = power_hi),
    width = 0.02,           # skinny caps
    linewidth = 0.5,
    position = pos
  ) +
  scale_x_log10(
    breaks = budgets,
    labels = c("16",  "128",  "512",  "4096", "32768")
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  ggplot2::labs(
    x = "N_train",
    y = "Power",
    color = "Method",
    shape = "Method"
  ) +
  # colorblind-friendly palette
  viridis::scale_color_viridis(discrete = TRUE, end = 0.9) +
  # clean theme: no grid, crisp axes
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.ticks.length = grid::unit(2.5, "pt"),
    plot.caption = ggplot2::element_text(hjust = 1, margin = ggplot2::margin(t = 6))
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, size = 3)),
    shape = ggplot2::guide_legend(override.aes = list(size = 3))
  )

print(p_power)

# make the plotting panel square
p_power = p_power + ggplot2::theme(aspect.ratio = 1)

# choose one size for both outputs
fig_size_in = 6  # inches

ggplot2::ggsave(
  "results/drpt_power_vs_budget.png",
  p_power,
  width = fig_size_in, height = fig_size_in, units = "in", dpi = 300
)

ggplot2::ggsave(
  "results/drpt_power_vs_budget.pdf",
  p_power,
  width = fig_size_in, height = fig_size_in, units = "in",
  device = cairo_pdf
)

cat("Saved:\n",
    "- raw replicate results to results/drpt_trials.csv\n",
    "- summary to results/drpt_summary.csv\n",
    "- plots to results/ (PNG + PDF)\n")

#--------------------------------------------------------------------
# Combine all figure in one
#--------------------------------------------------------------------

# Combine pictures
library(magick)
files = c("results/drpt_power_vs_budget.png", 
          "results/roc_curves_weighted.png")
imgs = image_read(files)

# Pick a target height (use the smaller one to avoid upscaling artifacts)
target_h = min(image_info(imgs)$height)

# Scale each image to that exact height, preserving aspect ratio
imgs = image_scale(imgs, paste0("x", target_h))

# Draw (a)/(b)/(c)
label_draw = function(img, lab, x = 15, y = 35) {
  info = image_info(img)
  h = info$height
  cex_val = max(1, h / 400) 
  
  img = image_draw(img)  
  par(mar = c(0,0,0,0), xpd = NA, family = "HersheySans")
  # soft white box behind text
  w  = strwidth(lab, cex = cex_val)
  ht = strheight(lab, cex = cex_val)
  rect(x - 6, y - 1.6*ht, x + 6 + 1.1*w, y + 6,
       col = rgb(1,1,1,0.7), border = NA)
  text(x, y, labels = lab, cex = cex_val, col = "black", adj = c(0, 1))
  dev.off()
  img  
}

a = label_draw(imgs[1], "(a)")
b = label_draw(imgs[2], "(b)")

# Combine horizontally
combo = image_append(image_join(a, b), stack = FALSE)
image_write(combo, "results/combined_3.png")