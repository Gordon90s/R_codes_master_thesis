# RPD VS EVT VS RTD_EVT for 3 heavy tailed bivariate distribution

print = F

# ===================================
# CAUCHY DISTRIBUTION

n <- 1000
mc <- 100
kappa.range <- 10:500
range.kappa.estimation <- 100:250 # range over which gamma.estimator.fct is averaged | needs to be within k.range
ylim <- c(0.45,1.45)
true.gamma <- 1

set.seed(1)
data <- data.polar <- array(NA, dim = c(n, 2, mc))
norm.data <- matrix(NA, n, mc)
for(j in 1:mc){
  data.polar[,,j] <- rcauchy(n)
  data[,,j] <- cbind(data.polar[,,j][,1]*cos(data.polar[,,j][,2]), data.polar[,,j][,1]*sin(data.polar[,,j][,2]))
  norm.data[,j] <- sqrt(rowSums(data[,,j]*data[,,j]))  
}

data.cauchy <- data[,,1]


# estimate gamma and kappa from mc number of monte carlo simulations for an heavy tailed distribution
# applied to normed data
estimate.gamma.mc <- function(norm.data, kappa.range, range.kappa.estimation, 
                              ylim = NULL, hill.plot = T, true.gamma = NULL, plot.legend = F, cex = 1){

  
gamma.estimator.fct <- matrix(NA, nrow = mc, ncol = max(kappa.range))
par(mfrow = c(3,3), mar = c(4.2, 4.2, 1.3, 2), cex = cex)

# calculate for all k projections with EVT modification
for(j in 1:mc){
  
  # taking absolute values as data is symmetric around zero
  # for other see for example using
  # apply(scalars,2,median)

  if(j<10){
    gamma.estimator.fct[j,kappa.range] <- gamma.estimator(norm.data[,j],k.range = kappa.range, legend = F, 
                                       legend.position = "bottomright", hill.plot = hill.plot, 
                                       main = paste0("m ",j), xlab = expression(kappa), 
                                       estimate.hill = T, estimate.moment = F, ylim = ylim, ylab = expression(widehat(gamma)["H"]*"("*kappa*")"))$gamma.hill
    abline(v = min(range.kappa.estimation), lty = 2, col = "red")
    abline(v = max(range.kappa.estimation), lty = 2, col = "red")
    if(plot.legend == T){
      if(j == 1) text(min(range.kappa.estimation)+50,0.6, label = expression(kappa["min"]), col = "red")
      if(j == 1) text(max(range.kappa.estimation)+50,0.6, label = expression(kappa["max"]), col = "red")
    }
    if(!is.null(true.gamma)){
    abline(h = true.gamma, lty = 2)
    if(j == 1 & plot.legend == T) text(160, mean(gamma.estimator.fct[j,range.kappa.estimation]) + 0.13, label = bquote(gamma*" = "*.(true.gamma)), col = "black")  
      }
    abline(h = mean(gamma.estimator.fct[j,range.kappa.estimation]), lty = 2, lwd = 2, col = "red")
    if(j == 1 & plot.legend == T) text(355, mean(gamma.estimator.fct[j,range.kappa.estimation]) - 0.1, label = expression(widehat(gamma)["H, mean"]), col = "red")
  }else{
  
  
    gamma.estimator.fct[j,kappa.range] <- gamma.estimator(norm.data[,j],k.range = kappa.range, hill.plot = F, 
                              estimate.hill = T, estimate.moment = F)$gamma.hill
  }
}
gamma.estimator.final <- mean(apply(gamma.estimator.fct[,range.kappa.estimation],1,mean))
kappa <- round(mean(range.kappa.estimation))

results <- list(gamma = gamma.estimator.final, kappa = kappa) 
return(results)
}

cex = 0.9
if(print == T){pdf("gamma_estimation_rcauchy.pdf", width = 9, height = 7.9)}#end print

cauchy.results <- estimate.gamma.mc(norm.data, kappa.range, range.kappa.estimation, 
                                    ylim = ylim, true.gamma = true.gamma, plot.legend = T, cex = cex)

if(print == T){dev.off()}

cauchy.results 

# ===================================
# CLOVER DISTRIBUTION

n <- 1000
mc <- 100
kappa.range <- 30:500
range.kappa.estimation <- 150:300 # range over which gamma.estimator.fct is averaged | needs to be within k.range
ylim <- c(0.22,0.43)
true.gamma <- 1/3

data.polar <- array(NA, dim = c(n, 2, mc))
data <- data.polar <- array(NA, dim = c(n, 2, mc))
norm.data <- matrix(NA, n, mc)
for(j in 1:mc){
  data.polar[,,j] <- rclover(n)
  data[,,j] <- cbind(data.polar[,,j][,1]*cos(data.polar[,,j][,2]), data.polar[,,j][,1]*sin(data.polar[,,j][,2]))
  norm.data[,j] <- sqrt(rowSums(data[,,j]*data[,,j]))  
}
data.clover <- data[,,j]

set.seed(1)

# redefine function only for plotting parameters... other than that, identical to above
estimate.gamma.mc <- function(norm.data, kappa.range, range.kappa.estimation, 
                              ylim = NULL, hill.plot = T, true.gamma = NULL, plot.legend = F, cex = 1){

  
gamma.estimator.fct <- matrix(NA, nrow = mc, ncol = max(kappa.range))
par(mfrow = c(3,3), mar = c(4.2, 4.2, 1.3, 2), cex = cex)

# calculate for all k projections with EVT modification
for(j in 1:mc){
  
  # taking absolute values as data is symmetric around zero
  # for other see for example using
  # apply(scalars,2,median)

  if(j<10){
    gamma.estimator.fct[j,kappa.range] <- gamma.estimator(norm.data[,j],k.range = kappa.range, legend = F, 
                                       legend.position = "bottomright", hill.plot = hill.plot, 
                                       main = paste0("m ",j), xlab = expression(kappa), 
                                       estimate.hill = T, estimate.moment = F, ylim = ylim, ylab = expression(widehat(gamma)["H"]*"("*kappa*")"))$gamma.hill
    abline(v = min(range.kappa.estimation), lty = 2, col = "red")
    abline(v = max(range.kappa.estimation), lty = 2, col = "red")
    if(plot.legend == T){
      if(j == 1) text(min(range.kappa.estimation)+50,0.25, label = expression(kappa["min"]), col = "red")
      if(j == 1) text(max(range.kappa.estimation)+50,0.25, label = expression(kappa["max"]), col = "red")
    }
    if(!is.null(true.gamma)){
    abline(h = true.gamma, lty = 2)
    if(j == 1 & plot.legend == T) text(min(range.kappa.estimation)+80, mean(gamma.estimator.fct[j,range.kappa.estimation]) + 0.05, label = bquote(gamma*" = 1/3"), col = "black")  
      }
    abline(h = mean(gamma.estimator.fct[j,range.kappa.estimation]), lty = 2, lwd = 2, col = "red")
    if(j == 1 & plot.legend == T) text(390, mean(gamma.estimator.fct[j,range.kappa.estimation]) - 0.025, label = expression(widehat(gamma)["H, mean"]), col = "red")
  }else{
  
  
    gamma.estimator.fct[j,kappa.range] <- gamma.estimator(norm.data[,j],k.range = kappa.range, hill.plot = F, 
                              estimate.hill = T, estimate.moment = F)$gamma.hill
  }
}
gamma.estimator.final <- mean(apply(gamma.estimator.fct[,range.kappa.estimation],1,mean))
kappa <- round(mean(range.kappa.estimation))

results <- list(gamma = gamma.estimator.final, kappa = kappa) 
return(results)
}

cex = 0.9
if(print == T){pdf("gamma_estimation_rclover.pdf", width = 9, height = 7.9)}#end print
clover.results <- estimate.gamma.mc(norm.data, kappa.range, range.kappa.estimation, ylim = ylim, true.gamma = true.gamma, cex = cex, plot.legend = T)

if(print == T){dev.off()}
clover.results




# ===================================
# ELLIPTICAL DISTRIBUTION DISTRIBUTION

n <- 1000
mc <- 100
kappa.range <- 30:500
range.kappa.estimation <- 150:300 # range over which gamma.estimator.fct is averaged | needs to be within k.range
ylim <- c(0.2,0.45)
true.gamma <- 1/3

set.seed(1)
data <- array(NA, dim = c(n, 2, mc))
norm.data <- matrix(NA, n, mc)
for(j in 1:mc){
  data[,,j] <- relliptical(n)
  norm.data[,j] <- sqrt(rowSums(data[,,j]*data[,,j]))  
}

data.elliptical <- data[,,1]

if(print == T){pdf("gamma_estimation_relliptical.pdf", width = 9, height = 7.9)}#end print
elliptica.results <- estimate.gamma.mc(norm.data, kappa.range, range.kappa.estimation, ylim = ylim, true.gamma = true.gamma, cex = cex)


if(print == T){dev.off()}
elliptica.results


# ===================================
# CAUCHY DISTRIBUTION

set.seed(1)
depth <- RTD.EVT_and_RPD(data.cauchy, data.cauchy, n.proj = 1000, k = cauchy.results$kappa, gamma = cauchy.results$gamma)

# plot(depth$RPD[depth.order.RPD][most.extremes], depth$RTD.EVT[depth.order.RPD][most.extremes], pch = ".")
# cor(depth$RPD[depth.order.RPD][most.extremes], depth$RTD.EVT[depth.order.RPD][most.extremes])

if(print == T){pdf("DD_plots_rcauchy.pdf", width = 9, height = 5.7)}

cex = 0.9
par(mfrow = c(2,3), mar = c(4.2, 4, 1, 2) + 0.1, cex = cex)
plot(depth$RPD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.1,0.4, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD.EVT),3))))
plot(depth$RPD, depth$RTD, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]))
text(0.1,0.4, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD),3))))
plot(depth$RTD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.1,0.4, labels = bquote(rho*" = "*.(round(cor(depth$RTD, depth$RTD.EVT),3))))

plot.lim <- c(0,0.025)
plot(depth$RPD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RPD, depth$RTD, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RTD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
if(print == T){dev.off()}

# depth.order.RPD <- order(depth$RPD)
# most.extremes <- 1:20
# plot(depth$RPD[depth.order.RPD][most.extremes], depth$RTD.EVT[depth.order.RPD][most.extremes], pch = "*")
# plot(depth$RPD[depth.order.RPD][most.extremes], depth$RTD[depth.order.RPD][most.extremes], pch = "*")
# plot(depth$RTD[depth.order.RPD][most.extremes], depth$RTD.EVT[depth.order.RPD][most.extremes], pch = "*")


# ===================================
# CLOVER DISTRIBUTION

set.seed(1)
depth <- RTD.EVT_and_RPD(data.clover, data.cauchy, n.proj = 1000, k = cauchy.results$kappa, gamma = cauchy.results$gamma)

if(print == T){pdf("DD_plots_rclover.pdf", width = 9, height = 5.7)}

cex = 0.9
par(mfrow = c(2,3), mar = c(4.2, 4, 1, 2) + 0.1, cex = cex)
plot(depth$RPD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD.EVT),3))))
plot(depth$RPD, depth$RTD, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD),3))))
plot(depth$RTD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RTD, depth$RTD.EVT),3))))

plot.lim <- c(0,0.15)
plot(depth$RPD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RPD, depth$RTD, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RTD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
if(print == T){dev.off()}


# ===================================
# ELLIPTICAL DISTRIBUTION DISTRIBUTION

set.seed(1)
depth <- RTD.EVT_and_RPD(data.elliptical, data.cauchy, n.proj = 1000, k = cauchy.results$kappa, gamma = cauchy.results$gamma)


if(print == T){pdf("DD_plots_rellipitcal.pdf", width = 9, height = 5.7)}

cex = 0.9
par(mfrow = c(2,3), mar = c(4.2, 4, 1, 2) + 0.1, cex = cex)
plot(depth$RPD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD.EVT),3))))
plot(depth$RPD, depth$RTD, pch = ".", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RPD, depth$RTD),3))))
plot(depth$RTD, depth$RTD.EVT, pch = ".", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]))
text(0.15,0.31, labels = bquote(rho*" = "*.(round(cor(depth$RTD, depth$RTD.EVT),3))))

plot.lim <- c(0,0.15)
plot(depth$RPD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RPD, depth$RTD, pch = "*", xlab = expression("D"["RP,1000"]), ylab = expression("D"["RT,1000"]), xlim = plot.lim, ylim = plot.lim)
plot(depth$RTD, depth$RTD.EVT, pch = "*", xlab = expression("D"["RT,1000"]), ylab = expression("D"["RT,EVT,1000"]), xlim = plot.lim, ylim = plot.lim)
if(print == T){dev.off()}
