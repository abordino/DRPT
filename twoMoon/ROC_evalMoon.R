rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/twoMoon")
#setwd("/storage/stats/strtng/SPT/twoMoon")

library(MASS)
library(latex2exp)
library(DRPT)
library(kerTests)
library(reticulate)
library(ggplot2)

set.seed(221198)

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

d = 4 
budgets = c(2^4, 2^7, 2^9, 2^12, 2^15) 
methods = c("bnre", "nre")

# training and test sizes for the ROC diagnostic
N_train = 50000      # per class for training (you set N = 200 above)
N_test  = 50000  # per class for evaluation

dir.create("results", showWarnings = FALSE, recursive = TRUE)

roc_points_list = list()
auc_summary_list = list()

for (B in budgets) {
  message("Budget: ", B)
  for (m in methods) {
    message("  Method: ", m)
    
    model_path = sprintf("artifacts/moon_rhat_%s_B%d.pt", m, B)
    if (!file.exists(model_path)) {
      warning(sprintf("Missing model file: %s — skipping.", model_path))
      next
    }
    r_hat = rhat_loader(model_path)
    
    # ---------- TRAIN SET ----------
    # positives = joint, label 1, weight 1
    thJ_tr = rsim_theta(N_train)
    xJ_tr  = rsim_x_given_theta(thJ_tr)
    Y_tr   = cbind(thJ_tr, xJ_tr)
    yJ_tr  = rep(1L, N_train)
    wJ_tr  = rep(1,  N_train)
    
    # negatives = product, label 0, weight = r_hat(theta,x)
    thP_tr = rsim_theta(N_train)
    xP_tr  = rsim_x_given_theta(rsim_theta(N_train))
    X_tr   = cbind(thP_tr, xP_tr)
    yP_tr  = rep(0L, N_train)
    
    # ratio weights for product samples
    r_tr = r_hat(X_tr[,1], X_tr[,2], X_tr[,3], X_tr[,4])
    r_tr[!is.finite(r_tr)] = NA_real_
    # small floor to avoid zero/NA weights in fitting
    r_tr = pmax(r_tr, 1e-12)
    
    # assemble training frame
    train_df = rbind(
      data.frame(theta1 = Y_tr[,1], theta2 = Y_tr[,2], x1 = Y_tr[,3], x2 = Y_tr[,4],
                 y = yJ_tr, w = wJ_tr),
      data.frame(theta1 = X_tr[,1], theta2 = X_tr[,2], x1 = X_tr[,3], x2 = X_tr[,4],
                 y = yP_tr, w = r_tr)
    )
    train_df = subset(train_df, is.finite(w) & !is.na(w))
    
    # Weighted logistic regression (scores are monotone with optimal test)
    clf = glm(y ~ theta1 + theta2 + x1 + x2,
               data = train_df,
               family = binomial(),
               weights = w)
    
    # ---------- TEST SET ----------
    # fresh positives (joint)
    thJ_te = rsim_theta(N_test)
    xJ_te  = rsim_x_given_theta(thJ_te)
    Y_te   = cbind(thJ_te, xJ_te)
    yJ_te  = rep(1L, N_test)
    wJ_te  = rep(1,  N_test)
    
    # fresh negatives (product) with weights = r_hat
    thP_te = rsim_theta(N_test)
    xP_te  = rsim_x_given_theta(rsim_theta(N_test))
    X_te   = cbind(thP_te, xP_te)
    yP_te  = rep(0L, N_test)
    r_te   = r_hat(X_te[,1], X_te[,2], X_te[,3], X_te[,4])
    r_te[!is.finite(r_te)] = NA_real_
    r_te   = pmax(r_te, 1e-12)
    
    test_df = rbind(
      data.frame(theta1 = Y_te[,1], theta2 = Y_te[,2], x1 = Y_te[,3], x2 = Y_te[,4],
                 y = yJ_te, w = wJ_te),
      data.frame(theta1 = X_te[,1], theta2 = X_te[,2], x1 = X_te[,3], x2 = X_te[,4],
                 y = yP_te, w = r_te)
    )
    test_df = subset(test_df, is.finite(w) & !is.na(w))
    
    # classifier scores on test (use logit/linear predictor; ROC is invariant to monotone transforms)
    score = as.numeric(predict(clf, newdata = test_df, type = "link"))
    
    # ---------- Weighted ROC (weighted by 1 for joint, r_hat for product) ----------
    tp_fp = WeightedROC::WeightedROC(
      guess  = score,
      label  = as.integer(test_df$y),  # 0/1
      weight = test_df$w
    )
    auc_val = as.numeric(WeightedROC::WeightedAUC(tp_fp))
    
    roc_points_list[[length(roc_points_list) + 1]] = transform(
      tp_fp,
      method = toupper(m),
      budget = B
    )
    auc_summary_list[[length(auc_summary_list) + 1]] = data.frame(
      method = toupper(m),
      budget = B,
      auc    = auc_val,
      n_eval = nrow(test_df)
    )
    message(sprintf("    AUC (weighted) = %.4f", auc_val))
  }
}

# ---------- Load & Plot (from saved files) ----------
suppressPackageStartupMessages({
  library(ggplot2)
  library(latex2exp)
})

# 1) Load
roc_points  = utils::read.csv("results/roc_points_weighted.csv",  stringsAsFactors = FALSE)
auc_summary = utils::read.csv("results/roc_auc_weighted.csv",    stringsAsFactors = FALSE)

# 2) Basic checks
needed = c("method","budget","FPR","TPR")
missing_cols = setdiff(needed, names(roc_points))
if (length(missing_cols)) stop("Missing columns in roc_points: ", paste(missing_cols, collapse=", "))

# 3) Coerce types (CSV sometimes loads numbers as character)
roc_points$method = as.character(roc_points$method)
roc_points$budget = suppressWarnings(as.numeric(roc_points$budget))
roc_points$FPR    = suppressWarnings(as.numeric(roc_points$FPR))
roc_points$TPR    = suppressWarnings(as.numeric(roc_points$TPR))

# 4) Budgets/methods
budgets = sort(unique(roc_points$budget))
methods_present = unique(roc_points$method)

# Nice budget factor (not strictly needed for plotting)
roc_points$budget_f = factor(roc_points$budget,
                              levels = budgets,
                              labels = formatC(budgets, format = "d", big.mark = ","))

# 5) Colors and legend levels (works if only one method exists)
bnre_base = "#0072B2"  # blue
nre_base  = "#D55E00"  # orange
bnre_grad = grDevices::colorRampPalette(c("#E6F0FA", bnre_base))(length(budgets))
nre_grad  = grDevices::colorRampPalette(c("#FBE6D8", nre_base))(length(budgets))
budget_lab = formatC(budgets, format = "d", big.mark = ",")

curve_levels = character(0)
curve_cols   = character(0)
if ("BNRE" %in% methods_present) {
  curve_levels = c(curve_levels, paste0("BNRE — ", budget_lab))
  curve_cols   = c(curve_cols,   bnre_grad)
}
if ("NRE" %in% methods_present) {
  curve_levels = c(curve_levels, paste0("NRE — ",  budget_lab))
  curve_cols   = c(curve_cols,   nre_grad)
}
# Fallback if methods are different labels
others = setdiff(methods_present, c("BNRE","NRE"))
if (length(others)) {
  # simple grayscale ramp per extra method
  for (m in others) {
    greys = grDevices::colorRampPalette(c("#F2F2F2", "#4D4D4D"))(length(budgets))
    curve_levels = c(curve_levels, paste0(m, " — ", budget_lab))
    curve_cols   = c(curve_cols,   greys)
  }
}
names(curve_cols) = curve_levels

# 6) Label per curve (method × budget) and lock factor order
roc_points$curve_label = with(
  roc_points,
  paste0(method, " — ", formatC(budget, format = "d", big.mark = ","))
)
roc_points$curve_label = factor(roc_points$curve_label, levels = curve_levels)

# 7) Plot
p_roc = ggplot(
  roc_points,
  aes(x = FPR, y = TPR,
      color = curve_label,
      group = interaction(method, budget))
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5,
              linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9, alpha = 0.95, na.rm = TRUE) +
  scale_color_manual(values = curve_cols, name = "Method — budget", drop = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title    = paste0("Weighted ROC: Joint vs Product reweighted by r̂"),
    subtitle = paste0("Classifier weights: w = 1 (joint), w = r̂",
                       " (product-of-marginals)"),
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10.5),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 1.2, alpha = 1)))

# 8) Save safely (create dir; force portable device)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
print(p_roc)
ggsave(filename = "results/roc_curves_weighted.png",
       plot = p_roc, width = 7.2, height = 5.2, dpi = 200)
