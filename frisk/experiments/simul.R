rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/distr_shift/simulationCpp/frisk")

library(latex2exp)
library(rootSolve)

source(file = "discreteDRPT.R")

###########################################################################
##### process the data 
###########################################################################

data.frisk11 = read.csv(file = "experiments/2015.csv")[, c(
  "pistol", "riflshot", "asltweap", "knifcuti",
  "machgun", "othrweap", "race", "crimsusp"
)]
data.frisk12 = read.csv(file = "experiments/2016.csv")[, c(
  "pistol", "riflshot", "asltweap", "knifcuti",
  "machgun", "othrweap", "race", "crimsusp"
)]
data = rbind(data.frisk11, data.frisk12)
data[1:2,]

# https://medium.com/@daniellekutner/stop-and-frisk-data-visualization-and-analysis-504c9a41ab6c
# https://github.com/VeeLeeKoh/Stop-and-Frisk-Data-Analysis

data = subset(data, race %in% c("W", "B"))
data = subset(data, crimsusp %in% c("CPW"))
data = na.omit(data)

# Map "Y" to 1 and "N" to 0
guilt_cols = c(
  "pistol", "riflshot", "asltweap", "knifcuti", 
  "machgun", "othrweap"
)

data[guilt_cols] = lapply(data[guilt_cols], function(col) {
  ifelse(col == "Y", 1, 0)
})

data$guilty = apply(data[guilt_cols], 1, max, na.rm = TRUE)
data = data[, !(names(data) %in% guilt_cols)]

# Subset for race == {P, Q (Hispanic) / B (Black)}
data_B = subset(data, race %in% c("B"))
# data_B = as.numeric(data_B[, !(names(data_B) %in% c("race"))])
data_B = as.numeric(data_B[, !(names(data_B) %in% c("race", "crimsusp"))])

# Subset for race == "W"
data_W = subset(data, race == "W")
# data_W = as.numeric(data_W[, !(names(data_W) %in% c("race"))])
data_W = as.numeric(data_W[, !(names(data_W) %in% c("race", "crimsusp"))])

###########################################################################
# Final dataset are data_W and data_B. Check ratio is nearly 5
###########################################################################

TB = table(data_B)
TW = table(data_W)

b1 = TB[2]/sum(TB)
w1 = TW[2]/sum(TW)
w1/b1

###########################################################################
## run discrete DRPT 
###########################################################################
xxx = 70
pval_theory = c()
pval_l2 = c()

interval = seq(from = 0.1, to = 7, length.out = xxx)
for (r2 in interval){
  print(r2)
  r = c(1, r2)
  r.func = function(x){
    return(r[x+1])
  }
  
  # compute hat.lambda for the test statistic
  Z = rbind(as.matrix(data_W),as.matrix(data_B))
  n = length(data_B); m = length(data_W)
  sum_lambda = function(l) {
    sum = 0
    for (k in 1:(n+m)) {
      Zk = Z[k,,drop = FALSE]
      sum = sum + 1 / (n + m * l * r.func(Zk))
    }
    return(sum - 1)
  }
  lambda.star = uniroot.all(sum_lambda, c(0, 100), tol = (.Machine$double.eps)^4)[1]
  
  
  pval_theory = c(pval_theory, discrete.DRPT(X = data_B, Y = data_W, r, H = 99) )
  pval_l2 = c(pval_l2, discrete.DRPT(X = data_B, Y = data_W, r, H = 99, 
                                     type = "V", lambda = lambda.star) )
}

###########################################################################
## Save data into .csv
###########################################################################
data_to_save = data.frame(Lambda = interval, 
                          pval_theory = pval_theory, pval_l2 = pval_l2)
write.csv(data_to_save, "experiments/results/frisk_pval15-16.csv", row.names = FALSE)


###########################################################################
## Plot results
###########################################################################
data12 = read.csv("experiments/results/frisk_pval11-12.csv")
data16 = read.csv("experiments/results/frisk_pval15-16.csv")
png("experiments/pictures/friskAppendix.png")
par(mfrow = c(1, 1))
plot(interval,
     data12$pval_l2, type = "b", col="blue", ylim = c(0,1),
     ylab='p-value', xlab=TeX(r'($r$)'), pch = 15, cex = 0.8)
points(interval,
     data16$pval_l2, type = "b", col="brown", ylim = c(0,1),
     ylab='p-value', xlab=TeX(r'($r$)'), pch = 15, cex = 0.8)
points(interval,
     data12$pval_theory, type = "b", col="cyan", ylim = c(0,1),
     ylab='p-value', xlab=TeX(r'($r$)'), pch = 15, cex = 0.8)
points(interval,
       data16$pval_theory, type = "b", col="pink", ylim = c(0,1),
       ylab='p-value', xlab=TeX(r'($r$)'), pch = 15, cex = 0.8)
abline(h = 0.05, col="red", lty = 2)
# legend("topleft", inset = 0.005,
#        horiz = FALSE, lty = 1, bty = "n",
#        legend = c("2011/2012 - V", "2015/2016 - V"),
#        col = c("blue", "brown"),
#        pch = c(15,15))
legend("topleft", inset = 0.001,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("2011/2012 - V", "2015/2016 - V", "2011/2012 - D", "2015/2016 - D"),
       col = c("blue", "brown", "cyan", "pink"),
       pch = c(15,15, 15, 15))
dev.off()
