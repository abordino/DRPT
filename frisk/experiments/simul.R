rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/frisk")

library(latex2exp)
library(rootSolve)
library(DRPT)

set.seed(120720)

# ----------------------------
# Datasets can be found at https://www.nyc.gov/site/nypd/stats/reports-analysis/stopfrisk.page
# ----------------------------

# ----------------------------
# Helpers to load/clean a pair
# ----------------------------
prep_years = function(file1, file2) {
  keep = c("pistol","riflshot","asltweap","knifcuti","machgun","othrweap","race","crimsusp")
  guilt_cols = c("pistol","riflshot","asltweap","knifcuti","machgun","othrweap")
  d1 = read.csv(file1)[, keep]
  d2 = read.csv(file2)[, keep]
  d  = rbind(d1, d2)
  d  = subset(d, race %in% c("W","B"))
  d  = subset(d, crimsusp %in% c("CPW"))
  d  = na.omit(d)
  for (col in guilt_cols) d[[col]] = ifelse(d[[col]] == "Y", 1, 0)
  d$guilty = apply(d[guilt_cols], 1, max, na.rm = TRUE)
  d = d[, !(names(d) %in% guilt_cols)]
  yB = as.numeric(subset(d, race == "B")[, !(names(d) %in% c("race","crimsusp"))])
  yW = as.numeric(subset(d, race == "W")[, !(names(d) %in% c("race","crimsusp"))])
  list(yW = yW, yB = yB)
}

# --------------------------------
# Wald OR + 95% CI (H-A correction)
# --------------------------------
wald_or_ci = function(yW, yB, alpha = 0.05, correction_if_zero = TRUE) {
  TW = table(factor(yW, levels = c(0,1)))
  TB = table(factor(yB, levels = c(0,1)))
  a = as.numeric(TW["1"])  # W ones
  b = as.numeric(TW["0"])  # W zeros
  c = as.numeric(TB["1"])  # B ones
  d = as.numeric(TB["0"])  # B zeros
  if (any(c(a,b,c,d) == 0) && correction_if_zero) {
    a = a + 0.5; b = b + 0.5; c = c + 0.5; d = d + 0.5   # Haldane–Anscombe
  }
  rhat   = (a*d)/(b*c)
  se_log = sqrt(1/a + 1/b + 1/c + 1/d)
  z      = qnorm(1 - alpha/2)
  ci     = exp(log(rhat) + c(-1, 1) * z * se_log)
  list(rhat = rhat, ci = ci, counts = c(a=a,b=b,c=c,d=d))
}

# --------------------------------------------
# Compute both DRPT p-value curves for a pair
# --------------------------------------------
compute_pair_curves = function(file_y1, file_y2, pair_tag,
                                xxx = 70, H = 99,
                                out_csv = sprintf("experiments/results/frisk_pvals_%s.csv", pair_tag)) {
  dat = prep_years(file_y1, file_y2)
  data_W = dat$yW
  data_B = dat$yB
  
  # Grid for r = (1, r2)
  interval = seq(from = 0.1, to = 7, length.out = xxx)
  
  # Precompute bits used in lambda root for type="V"
  Z = rbind(as.matrix(data_W), as.matrix(data_B))
  n = length(data_B); m = length(data_W)
  
  r.func = function(x, r) r[x + 1]
  
  pval_V = numeric(length(interval))
  pval_D = numeric(length(interval))
  
  for (i in seq_along(interval)) {
    r2 = interval[i]
    r  = c(1, r2)
    
    # ---- DRPT p-values ----
    pD = tryCatch(
      discrete.DRPT(X = data_B, Y = data_W, r = r, H = H, type = "D"),
      error = function(e) NA_real_
    )
    pV = tryCatch(
      discrete.DRPT(X = data_B, Y = data_W, r = r, H = H, type = "V"),
      error = function(e) NA_real_
    )
    
    pval_D[i] = pD
    pval_V[i] = pV
  }
  
  dir.create("experiments/results", showWarnings = FALSE, recursive = TRUE)
  out_df = data.frame(Lambda = interval, pval_D = pval_D, pval_V = pval_V)
  write.csv(out_df, out_csv, row.names = FALSE)
  
  # return also the cleaned data for Wald
  list(curves = out_df, yW = data_W, yB = data_B)
}

# -------------------------------------------------
# Compute for both pairs + Wald, then 4-curve plot
# -------------------------------------------------
res_1112 = compute_pair_curves("experiments/2011.csv", "experiments/2012.csv", "11-12")
res_1516 = compute_pair_curves("experiments/2015.csv", "experiments/2016.csv", "15-16")

data_1112 = res_1112$curves
data_1516 = res_1516$curves

wald_1112 = wald_or_ci(res_1112$yW, res_1112$yB)
wald_1516 = wald_or_ci(res_1516$yW, res_1516$yB)

# Common x grid (identical by construction)
xgrid = data_1112$Lambda

# x-limits to include CI endpoints and point estimates
xlim_rng = range(c(xgrid, wald_1112$ci, wald_1516$ci, wald_1112$rhat, wald_1516$rhat), finite = TRUE)

# ----------------------
# Plot: 4 curves + Wald
# ----------------------
dir.create("experiments/pictures", showWarnings = FALSE, recursive = TRUE)
png("experiments/pictures/frisk_4curves_wald.png")

# Consistent styling for both curves
curve_lwd = 2.5
curve_lty = 1     # solid
pch_1112  = 16    # filled circle
pch_1516  = 17    # filled triangle

plot(xgrid, data_1112$pval_V,
     type = "b", lwd = curve_lwd, lty = curve_lty, pch = pch_1112,
     col = "blue", ylim = c(0, 1), xlim = xlim_rng,
     xlab = TeX(r'($r$)'), ylab = "p-value",
     xaxs = "i", yaxs = "i")

lines(xgrid, data_1516$pval_V,
      type = "b", lwd = curve_lwd, lty = curve_lty, pch = pch_1516,
      col = "brown")

abline(h = 0.05, col = "red", lty = 2)  # dashed threshold

# Wald overlays
abline(v = wald_1112$rhat, col = "blue",  lty = 3)   # dotted verticals
abline(v = wald_1516$rhat, col = "brown", lty = 3)

# CI as horizontal segments (slightly staggered y to avoid overlap)
y_ci_1112 = 0.50
y_ci_1516 = 0.55
segments(wald_1112$ci[1], y_ci_1112, wald_1112$ci[2], y_ci_1112, col = "blue",  lwd = 2)
segments(wald_1516$ci[1], y_ci_1516, wald_1516$ci[2], y_ci_1516, col = "brown", lwd = 2)
points(wald_1112$rhat, y_ci_1112, pch = 16, col = "blue",  cex = 0.9)
points(wald_1516$rhat, y_ci_1516, pch = 16, col = "brown", cex = 0.9)

legend("topright", inset = c(0.21, 0), bty = "n",
       legend = c("2011/2012 – (E3)", "2015/2016 – (E3)",
                  "Wald mean (vertical)", "Wald 95% CI"),
       col = c("blue","brown","black","black"),
       lwd = c(curve_lwd, curve_lwd, 1, 2),
       lty = c(curve_lty, curve_lty, 3, 1),
       pch = c(pch_1112, pch_1516, NA, NA))

dev.off()


---------------------------------------------------------------
  
  # Combine pictures
  library(magick)
  files = c("experiments/pictures/binaryRvarying copy.png",
            "experiments/pictures/4d_combinedAppendix copy.png" 
          )
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
image_write(combo, "experiments/pictures/combined_4.png")