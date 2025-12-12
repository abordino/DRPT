rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/synthetic")

library(MASS)
library(latex2exp)
library(squash)
library(densratio)
library(mvtnorm)

source(file = "SPT.R")
source(file = "estimate_marginal_ratio3.R")

invCDF = function(x){
  return(sqrt(x))
}

####################################
############ uLSIF d = 1
####################################

set.seed(123)
n_train <- 200

# Generate samples d = 1
X <- runif(n_train)
Y <- invCDF(runif(n_train))

# Run KLIEP estimation
densratio_obj <- densratio(Y, X)

new_x <- seq(0, 1, by = 0.02)
estimated_density_ratio <- densratio_obj$compute_density_ratio(new_x)

plot(new_x, estimated_density_ratio, pch=19)
lines(new_x, 2*new_x, type="l", col = "red")


####################################
############ uLSIF d = 2
####################################
n_train <- 200

X <- data.frame(x1 = runif(n_train), x2 = runif(n_train))
Y <- data.frame(x1 = invCDF(runif(n_train)), x2 = runif(n_train))

densratio_obj_d2 <- densratio(Y, X)

true_density_ratio <- function(x) {
  return(2*x$Var1)
}

N <- 20
range <- seq(0, 1, length.out = N)
input <- expand.grid(range, range)
w_true <- matrix(true_density_ratio(input), nrow = N)
w_hat <- matrix(densratio_obj_d2$compute_density_ratio(input), nrow = N)

par(mfrow = c(1, 2))
contour(range, range, w_true, main = "True Density Ratio")
contour(range, range, w_hat, main = "Estimated Density Ratio")

###################################
## LL d = 1
###################################
n_train <- 20000

# Generate samples d = 1
X <- runif(n_train)
Y <- invCDF(runif(n_train))

n_x <- length(X)
n_y <- length(Y)

df_train <- data.frame(
  x      = c(X, Y),                          # predictor
  label  = factor(c(rep(0, n_x), rep(1, n_y)))  # 0 = Y, 1 = X
)

## logistic regression (linear log-odds model)
fit.marginal <- glm(label ~ x, data = df_train, family = binomial())

## function that returns \hat p(x)/\hat q(x)
marginal_ratio <- function(z){
  prob <- predict(fit.marginal,
                  newdata = data.frame(x = z),
                  type = "response")          # \hat P(label = 1 | x = z)
  prob / (1 - prob)                          # density ratio estimate
}


new_x <- seq(0, 1, by = 0.02)
estimated_density_ratio <- marginal_ratio(new_x)

plot(new_x, estimated_density_ratio, pch=19)
lines(new_x, 2*new_x, type="l", col = "red")