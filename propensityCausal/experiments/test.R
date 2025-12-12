# ==========================================
# p = 1 sanity: eta(x) = 0.05*x + gamma*sin(omega*x), omega=3
# Truth: r(x) = exp(eta(x)),  pi = 0.5  (odd symmetry)
# Estimators: (i) LL -> ratio = exp(eta_hat), (ii) KLR (CVST + rbfdot)
# Saves PNG + CSV per gamma
# ==========================================

rm(list = ls()); gc()
set.seed(170966)

## ------------ Config ------------
gammas   <- c(0, 0.25, 0.5, 1, 2)
omega    <- 3                  # frequency in sin(omega * x)
N_train  <- 6000
x_grid   <- seq(-3.5, 3.5, length.out = 401)

# KLR kernel & regularization
# We'll pick sigma from a small data-driven grid; lambda fixed (can tune if needed)
klr_lambda <- 5e-4
sigma_multipliers <- c(0.5, 1, 2, 4)   # grid around robust scale

# Output
out_dir_plots <- "experiments/pictures"
out_dir_csv   <- "experiments/results"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_csv,   recursive = TRUE, showWarnings = FALSE)

## ------------ Deps ------------
pkgs <- c("MASS","CVST","kernlab")
need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## ------------ DGP (p = 1) ------------
expit   <- function(x) 1/(1+exp(-x))
eta_fun <- function(x, g) 0.05*x + g*sin(omega*x)    # odd in x -> pi = 0.5
draw_joint <- function(N, g) {
  x <- as.numeric(MASS::mvrnorm(N, mu = 0, Sigma = matrix(1,1,1)))
  e <- expit(eta_fun(x, g))
  D <- rbinom(N, 1, e)
  list(Z = matrix(x, ncol = 1), D = D)
}

## ------------ Ground truth (analytic) ------------
true_ratio <- function(g, grid = x_grid) {
  list(x = grid, r_true = exp(eta_fun(grid, g)), pi = 0.5)
}

## ------------ Estimator (i): Logistic (LL) ------------
# IMPORTANT: ratio(x) = exp(eta_hat(x) - logit(pi_target)); here pi_target=0.5 -> exp(eta_hat)
fit_ratio_LL <- function(N_train, g, pi_target = 0.5) {
  S <- draw_joint(N_train, g)
  df <- data.frame(D = as.numeric(S$D), z1 = as.numeric(S$Z[,1]))
  fit <- glm(D ~ z1, data = df, family = binomial())
  function(xnew) {
    eta_hat <- as.numeric(predict(fit, newdata = data.frame(z1 = as.numeric(xnew)), type = "link"))
    exp(eta_hat - qlogis(pi_target))  # with pi_target=0.5, this is just exp(eta_hat)
  }
}

## ------------ KLR helpers ------------
robust_scale_1d <- function(x) {
  # robust sd ~ MAD * 1.4826; fallback to sd; fallback to 1
  s <- median(abs(x - median(x))) * 1.4826
  if (!is.finite(s) || s <= 1e-12) s <- sd(x)
  if (!is.finite(s) || s <= 1e-12) s <- 1
  s
}

## ------------ Estimator (ii): KLR (CVST + rbfdot) ------------
# We select sigma from a small grid around the robust scale (data-driven),
# train KLR for each sigma, pick the one with highest in-sample log-likelihood (simple & fast).
fit_ratio_KLR <- function(N_train, g,
                          lambda = klr_lambda,
                          sigma_mult = sigma_multipliers,
                          tol = 1e-6, maxiter = 500) {
  S <- draw_joint(N_train, g)
  Z <- S$Z; D <- S$D
  x <- as.numeric(Z[,1])
  n0 <- sum(D==0); n1 <- sum(D==1)
  if (n0 < 10 || n1 < 10) stop("Too few samples per class for KLR.")
  
  # CVST data container (pooled)
  x.fit      <- matrix(x, ncol = 1)
  label.fit  <- factor(D)  # {0,1}
  data.fit   <- CVST::constructData(x.fit, label.fit)
  klrlearner <- CVST::constructKlogRegLearner()
  
  # Build sigma grid (rbf kernel uses parameter 'sigma' in kernlab::rbfdot)
  base <- robust_scale_1d(x)
  sigma_grid <- pmax(1e-6, base * sigma_mult)
  
  best <- NULL
  best_ll <- -Inf
  for (sg in sigma_grid) {
    params <- list(kernel = "rbfdot",
                   sigma  = sg,
                   lambda = lambda,
                   tol    = tol,
                   maxiter= maxiter)
    fit <- klrlearner$learn(data.fit, params)
    
    # compute in-sample log-likelihood quickly
    fx <- kernlab::kernelMult(fit$kernel, x.fit, fit$data, fit$alpha)
    p  <- 1/(1+exp(-as.vector(fx)))
    p  <- pmin(pmax(p, 1e-12), 1-1e-12)
    ll <- sum(D*log(p) + (1-D)*log(1-p))
    if (ll > best_ll) { best_ll <- ll; best <- fit }
  }
  
  # Return ratio function: r(x) = (n0/n1) * odds_KLR(x)
  function(xnew) {
    xnew <- matrix(as.numeric(xnew), ncol = 1)
    fx <- kernlab::kernelMult(best$kernel, xnew, best$data, best$alpha)
    p  <- 1/(1+exp(-as.vector(fx)))
    p  <- pmin(pmax(p, 1e-12), 1-1e-12)
    odds <- p/(1-p)
    as.numeric( (n0 / n1) * odds )
  }
}

## ------------ Run one gamma ------------
run_one_gamma <- function(g) {
  message(sprintf("==> gamma = %g", g))
  
  T1   <- true_ratio(g, x_grid)
  r_LL <- fit_ratio_LL(N_train, g, pi_target = 0.5)
  r_KL <- fit_ratio_KLR(N_train, g)
  
  y_LL <- r_LL(T1$x)
  y_KL <- r_KL(T1$x)
  
  # Save & plot
  df <- data.frame(x1 = T1$x, r_true = T1$r_true, r_LL = y_LL, r_KLR = y_KL)
  csv_path <- file.path(out_dir_csv, sprintf("sanity_ratio_p1_KLR_gamma_%g.csv", g))
  write.csv(df, csv_path, row.names = FALSE)
  
  png_path <- file.path(out_dir_plots, sprintf("sanity_ratio_p1_KLR_gamma_%g.png", g))
  png(png_path, width = 1000, height = 650)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  plot(df$x1, df$r_true, type = "l", lwd = 3,
       xlab = expression(x[1]), ylab = "density ratio  r(x1)",
       main = bquote(eta(x) == 0.05*x + .(g) * sin(.(omega)*x)))
  lines(df$x1, df$r_LL,  lwd = 2, lty = 2, col = "blue")   # LL
  lines(df$x1, df$r_KLR, lwd = 2, lty = 3, col = "red")    # KLR
  abline(h = 1, col = "gray70", lty = 3)
  legend("topleft", bty = "n", lwd = c(3,2,2), lty = c(1,2,3),
         col = c("black","blue","red"),
         legend = c("TRUE (analytic)", "LL", "KLR"))
  dev.off()
  
  message("  CSV:  ", csv_path)
  message("  Plot: ", png_path)
  invisible(df)
}

## ------------ Go ------------
invisible(lapply(gammas, run_one_gamma))
message("Done.")

