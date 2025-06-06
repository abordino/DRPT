rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/stroop")

library(MASS)
library(latex2exp)

source(file = "DRPT.R")


###########################################################################
## process the data 
###########################################################################

X = c(); Y = c()
data = read.csv("Stro.csv", header = TRUE)
for (i in 1:dim(data)[1]){
  str = strsplit(data[i,], " ")[[1]]
  X = c(X, as.numeric(str[2]))
  Y = c(Y, as.numeric(str[3]))
}

normalize = function(x) {
  (x - min(x)) / (max(x) - min(x))
}
X  = normalize(X); Y = normalize(Y)

d = 1
gaussian.kernel = function(x, y, lambda = 1){
  return(lambda^(-d) * exp(-sum(((x - y) ^ 2) / (lambda ^ 2))))
}

###########################################################################
## run DRPT 
###########################################################################

xxx = 30
pval_l = c()
for (l in seq(from = 0, to = 0.3, length.out = xxx)){
  print(l)
  
  r = function(x){
    return(exp(x/l))
  }
  
  pval_l = c(pval_l, DRPT(X,Y, S=50, H=99, r, gaussian.kernel))
}

###########################################################################
## Save data into .csv
###########################################################################
data_to_save = data.frame(Lambda = seq(from = 0, to = 0.3, length.out = xxx), 
                          P_lambda = pval_l)
write.csv(data_to_save, "experiments/results/stroop_pval.csv", row.names = FALSE)


###########################################################################
## Plot results
###########################################################################
data = read.csv("experiments/results/stroop_pval.csv")
png("experiments/pictures/p_val.png")
par(mfrow = c(1, 1))
plot(seq(from = 0, to = 0.3, length.out = xxx),
     data$P_lambda, type = "b", col="purple",
     ylab='p-value', xlab=TeX(r'($\eta$)'), pch = 14, cex = 0.8)
abline(h = 0.05, col="red", lty = 2)
legend("topleft", inset = 0.05,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("DRPT"),
       col = c("darkviolet"),
       pch = c(14))
dev.off()
