# ==========================================
# Causal DGP, LL ratio, DRPT (E1), save decisions
# ==========================================
rm(list = ls()); gc()

setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/propensityCausal/")
# setwd("/storage/stats/strtng/SPT/syntheticCausal/")

# --- deps ---
pkgs = c("MASS","DRPT","kerTests", "latex2exp")
need = pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(need) > 0) install.packages(need)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(170966)

# -------------------
# DGP definitions
# -------------------
p = 10
rho = 0.5
Sigma = outer(1:p, 1:p, function(i,j) rho^abs(i-j))  # AR(1)
beta0 = 0
beta  = c(0.6, -0.5, 0.4, 0.3, -0.2, 0.1, 0.25, -0.15, 0.2, -0.25)

expit = function(x) 1/(1+exp(-x))

eta_gamma = function(Z, gamma) {
  drop(beta0 + Z %*% beta) + gamma * ( sin(10 * Z[,1]) + Z[,2] * Z[,3] )
}

draw_joint = function(N, gamma) {
  Z = MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma)
  e = expit(eta_gamma(Z, gamma))
  D = rbinom(N, 1, e)
  list(Z = Z, D = D)
}

draw_groups_fixed = function(n0, n1, gamma) {
  Z0 = NULL; Z1 = NULL
  while (is.null(Z0) || nrow(Z0) < n0 || is.null(Z1) || nrow(Z1) < n1) {
    S = draw_joint(4000, gamma)
    Z0 = rbind(Z0, S$Z[S$D==0, , drop = FALSE])
    Z1 = rbind(Z1, S$Z[S$D==1, , drop = FALSE])
  }
  list(X = Z0[1:n0, , drop = FALSE],  # f = Z|D=0
       Y = Z1[1:n1, , drop = FALSE])  # g = Z|D=1
}

# -------------------
# Kernel for DRPT
# -------------------
gaussian.kernel = function(x, y, lambda = 1) {
  d = length(x)
  lambda^(-d) * exp(-sum((x - y)^2) / (lambda^2))
}

# -------------------
# Ratio estimator: Linear Logistic
# -------------------
fit_ratio_LL = function(N_train, gamma) {
  S = draw_joint(N_train, gamma)
  Z = S$Z; D = S$D
  df = as.data.frame(Z); colnames(df) = paste0("z", 1:ncol(Z))
  df$D = factor(D)
  pi_train = mean(D)
  
  fml = as.formula(paste("D ~", paste(colnames(df)[1:p], collapse = " + ")))
  glm_fit = glm(fml, data = df, family = binomial())
  
  # rhat(z) = exp(eta_hat(z) - logit(pi_train)) = ((1-pi)/pi)*odds
  rhat = function(...) {
    Znew = cbind(...)
    if (is.null(dim(Znew))) Znew = matrix(Znew, nrow = 1)
    Zdf = as.data.frame(Znew)
    colnames(Zdf) = paste0("z", 1:ncol(Znew))
    phat = predict(glm_fit, newdata = Zdf, type = "response")
    phat = pmin(pmax(phat, 1e-12), 1 - 1e-12)
    odds = phat / (1 - phat)
    exp(-qlogis(pi_train)) * odds
  }
  
  rhat
}

create_lookup_list = function(X, rX, Y, rY) {
  keysX = apply(X, 1, paste, collapse = "_")
  keysY = apply(Y, 1, paste, collapse = "_")
  
  lookup_list = setNames(as.list(c(rX, rY)), c(keysX, keysY))
  
  return(lookup_list)
}

# -------------------
# Experiment config
# -------------------
gammas = c(0, 0.25, 0.5, 1, 2)
N_train = 1000                 # training size for ratio
N_test_per_group = 150         # if n=m=100 => total 200 in test
MC = 30                       # repetitions per gamma
alpha_level = 0.05

# output dir
out_dir = file.path(getwd(), "experiments", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------
# Run experiment
# -------------------
results = list()
rowid = 1

for (g in gammas) {
  message(sprintf("Gamma = %s", g))
  
  # Fit ratio once per gamma on independent training data
  rfun = fit_ratio_LL(N_train = N_train, gamma = g)
  
  rej_count = 0
  
  for (b in seq_len(MC)) {
    G = draw_groups_fixed(n0 = N_test_per_group, n1 = N_test_per_group, gamma = g)
    X = G$X; Y = G$Y
    
    lookup_table = create_lookup_list(X, rfun(X), Y, rfun(Y))
    
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
    
    lam = kerTests::med_sigma(X, Y)
    kfun = function(x, y) gaussian.kernel(x, y, lambda = lam)
    
    pval = DRPT::DRPT(X, Y, r = r.hat, kernel = kfun)
    decision = as.integer(pval < alpha_level)
    
    results[[rowid]] = data.frame(gamma = g, rep = b,
                                   p_value = as.numeric(pval),
                                   decision = decision)
    rowid = rowid + 1
    rej_count = rej_count + decision
  }
  
  power_g = rej_count / MC
  message(sprintf("  -> Power estimate: %.3f", power_g))
}

res = do.call(rbind, results)

# ---- Save single-run decisions (and p-values)
runs_file = file.path(out_dir, "ll_drpt_runs.csv")
write.csv(res, runs_file, row.names = FALSE)

# Also save per-gamma decision sequences (rep,decision) if you like
for (g in gammas) {
  sub = res[res$gamma == g, c("rep","decision")]
  fn = file.path(out_dir, sprintf("ll_drpt_decisions_gamma_%s.csv", g))
  write.csv(sub, fn, row.names = FALSE)
}

# ---- Save power summary
power_tab = aggregate(decision ~ gamma, data = res, mean)
colnames(power_tab) = c("gamma", "power")
power_file = file.path(out_dir, "ll_drpt_power.csv")
write.csv(power_tab, power_file, row.names = FALSE)

print(power_tab)
message("Saved runs to: ", runs_file)
message("Saved power to: ", power_file)





# --- Plot (unchanged inputs) ---
data.LL  = read.csv(file = "experiments/results/ll_drpt_power.csv")
data.KLR   = read.csv(file = "experiments/results/klr_irls_drpt_power.csv")

# MC must be defined in your environment
MC = 100

compute_sd = function(power) {
  m = power
  sqrt(m * (1 - m) / MC)
}

sd_LL = compute_sd(data.LL$power)
sd_KLR  = compute_sd(data.KLR$power)

# Helper to draw error bars, clipped to [0,1]
add_err = function(x, y, sd, col) {
  ylow  = pmax(0, y - sd)
  yhigh = pmin(1, y + sd)
  arrows(x0 = x, y0 = ylow, x1 = x, y1 = yhigh,
         angle = 90, code = 3, length = 0.04,
         col = col, lwd = 1.3)
}

# Compute y-limits considering error bars
all_low  = pmin(data.LL$power  - sd_LL,
                data.KLR$power   - sd_KLR, na.rm = TRUE)
all_high = pmax(data.LL$power  + sd_LL,
                data.KLR$power   + sd_KLR, na.rm = TRUE)
ylim_use = c(max(0, min(all_low, na.rm = TRUE)),
             min(1, max(all_high, na.rm = TRUE)))

dir.create("experiments/pictures", recursive = TRUE, showWarnings = FALSE)
png("experiments/pictures/Propensity.png")
par(mfrow = c(1, 1))

plot(data.LL$gamma, data.LL$power, type = "b",
     col = "darkviolet", lwd = 2, pch = 14,
     xlab = TeX(r'($\gamma$)', bold = TRUE), ylab = "Power",
     ylim = ylim_use)

# Other series
lines(data.KLR$gamma, data.KLR$power,   
      type = "b", col = "orange", lwd = 2, pch = 16)

# SD segment bars
add_err(data.LL$gamma, data.LL$power,        sd_LL, "darkviolet")
add_err(data.KLR$gamma, data.KLR$power,        sd_KLR, "orange")

legend("topleft", inset = 0.01, horiz = FALSE, lty = 1, bty = "n",
       legend = c("(E1) - LL", "(E1) - KLR"),
       col = c("darkviolet", "orange"),
       pch = c(14, 16))

abline(h = 0.05, col = "red", lty = 2)
dev.off()