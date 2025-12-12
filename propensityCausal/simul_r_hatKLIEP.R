# ==========================================
# Causal DGP, KLR ratio (IRLS only, NO SCALING), DRPT (E1), save decisions
# ==========================================
rm(list = ls()); gc()

# --- (optional) working dir ---
setwd("~/Documents/phd/projects/DRPT/code/simulationCpp/propensityCausal/")
# setwd("/storage/stats/strtng/SPT/syntheticCausal/")

# --- deps ---
pkgs = c("MASS","DRPT","kerTests", "CVST", "densratio")
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

# exact group sizes for DRPT comparisons
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
# Kernel for DRPT (uses length-scale lambda)
# -------------------
gaussian.kernel = function(x, y, lambda = 1) {
  d = length(x)
  lambda^(-d) * exp(-sum((x - y)^2) / (lambda^2))
}

# --- Nonparametric ratio via uLSIF (g/f), cross-validated bandwidth ---
fit_ratio_KLIEP = function(N_train, gamma,
                           params = list(method = "KLIEP",
                                         kernel_num = 300,
                                         sigma = "auto",
                                         lambda = "auto")) {
  
  # draw and split
  S  = draw_joint(N_train, gamma)
  Z  = S$Z; D = S$D
  Xf = Z[D == 0, , drop = FALSE]   # f = Z|D=0  (denominator)
  Xg = Z[D == 1, , drop = FALSE]   # g = Z|D=1  (numerator)
  nf = nrow(Xf); ng = nrow(Xg)
  if (nf < 10 || ng < 10) return(NULL)
  
  # fit r(x) ≈ g(x)/f(x) with uLSIF
  dr = densratio::densratio(
    x1 = Xg, x2 = Xf,
    method = params$method,                   # "uLSIF" (or "KLIEP")
    kernel_num = min(params$kernel_num, ng),  # keep runtime stable
    sigma = params$sigma, lambda = params$lambda,
    verbose = FALSE
  )
  
  # callable r̂(x)
  rhat = function(...) {
    newX = cbind(...)
    if (is.null(dim(newX))) newX = matrix(newX, nrow = 1)
    out = dr$compute_density_ratio(newX)      # matrix OK
    pmin(pmax(as.numeric(out), 1e-12), 1e12)  # gentle clipping like before
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
N_test_per_group = 150         # n=m=100
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
  
  # Fit KLR ratio once per gamma on independent training data (IRLS, no scaling)
  rfun = fit_ratio_KLIEP(N_train = N_train, gamma = g)
  
  rej_count = 0
  valid_reps = 0
  
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
    
    # DRPT kernel bandwidth from RAW test features
    lam = kerTests::med_sigma(X, Y)
    kfun = function(x, y) gaussian.kernel(x, y, lambda = lam)
    
    pval = DRPT::DRPT(X, Y, r = rhat, kernel = kfun)
    decision = as.integer(pval < alpha_level)
    
    results[[rowid]] = data.frame(gamma = g, rep = b,
                                  p_value = as.numeric(pval),
                                  decision = decision)
    rowid = rowid + 1
  }
  
  power_g = if (valid_reps > 0) rej_count / valid_reps else NA_real_
  message(sprintf("  -> Power estimate (na.rm): %.3f over %d valid reps",
                  ifelse(is.na(power_g), NaN, power_g), valid_reps))
}

res = do.call(rbind, results)

# ---- Save single-run decisions (and p-values)
runs_file = file.path(out_dir, "kliep_drpt_runs.csv")
write.csv(res, runs_file, row.names = FALSE)

# Also save per-gamma decision sequences (rep,decision)
for (g in gammas) {
  sub = res[res$gamma == g, c("rep","decision")]
  fn = file.path(out_dir, sprintf("kliep_drpt_decisions_gamma_%s.csv", g))
  write.csv(sub, fn, row.names = FALSE)
}

# ---- Save power summary (ignore NA decisions)
ok = res[!is.na(res$decision), , drop = FALSE]
if (nrow(ok) > 0) {
  power_tab = aggregate(decision ~ gamma, data = ok, mean)
  colnames(power_tab) = c("gamma", "power")
  power_file = file.path(out_dir, "kliep_drpt_power.csv")
  write.csv(power_tab, power_file, row.names = FALSE)
  print(power_tab)
  message("Saved runs to: ", runs_file)
  message("Saved power to: ", power_file)
} else {
  message("No valid reps to summarize (all fits failed).")
}
