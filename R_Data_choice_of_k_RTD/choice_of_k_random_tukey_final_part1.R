# Choice of k RT

reload.packages = F

if(reload.packages == T){

#install.packages("ddalpha")
library(ddalpha)
#install.packages("clusterGeneration")
library(clusterGeneration)
#install.packages("xtable")
library(xtable) # to create LaTeX tables
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
}


# WARNING: work not (exactly) reproducible as seed function within "depth.halfspace" is not working. (as of 12.04.2017)


dim.mv.simulations.all <- c(2,3)
n.x.simulations.all <- c(50,100,250)

# n.proj.k.star <- 10
n.proj.k.star <- 10000



all.simulations <-  matrix(NA,length(dim.mv.simulations.all)*length(n.x.simulations.all),2)
for(i in 1:length(dim.mv.simulations.all)){
  for(j in 1:length(n.x.simulations.all)){
    all.simulations[j+(i-1)*length(n.x.simulations.all),] <- c(dim.mv.simulations.all[i],n.x.simulations.all[j])
    }
  }



# n.proj <- 10
n.proj <- 800


n.proj.seq <- 1:n.proj

result.matrix.final <- array(NA, dim = c(3,n.proj,3,nrow(all.simulations)))
k.zero.final <- matrix(NA, 3, nrow(all.simulations))

ptm <- proc.time() # measure time of the code

for(z2 in 1:nrow(all.simulations)){
print(z2)
ptm.z2 <- proc.time()
# smaller test

# number of MC simulations
mc <- 1000


# simulations
n.x.simulations <- all.simulations[z2,2] # slightly change in the algorithm, keep length of 1. (only for naming though, so no big deal)
dim.mv.simulations <- all.simulations[z2,1]

amount.of.simulations <- length(n.x.simulations)*length(dim.mv.simulations)
data.correlation.MC <- array(NA, dim = c(3,n.proj,amount.of.simulations,mc))
dimnames(data.correlation.MC)[[1]] <- c("k","RT","MH")

for(z in 1:mc){
# print(z)
if(1==1){


simulation.stats <- matrix(NA, nrow = amount.of.simulations, ncol = 3) # to save n.x, dim.mv and n.proj for each simulation
data.correlation <- array(NA, dim = c(3,n.proj,amount.of.simulations)) # to save the seq.sequences and the correlation values
dimnames(data.correlation)[[1]] <- c("k","RT","MH")


# set.seed <- 1 # Setting seeds does not anker the "depth.halfspace" function, so useless

for(l in 1:length(n.x.simulations)){
  for(m in 1:length(dim.mv.simulations)){

n.x <- n.x.simulations[l]
dim.mv <- dim.mv.simulations[m]


simulation.stats[(l-1)*length(dim.mv.simulations)+m, ] <- c(n.x, dim.mv, n.proj)

# print(paste("l =",l,", m = ",m,", n.sim = ",(l-1)*length(dim.mv.simulations)+m))
# print(paste("n.x =",n.x, "dim.mv =",dim.mv))


# simulate data (multivariate standard normal data)

cov.x <- diag(dim.mv)
x <- rmvnorm(n.x, mean = rep(0,dim.mv), sigma = cov.x)

# true Mahalanobis depth
depth.x.M.true <- 1/(1+mahalanobis(x, center = rep(0, ncol(x)), cov = cov.x))

#ptm <- proc.time() # measure time of the code

# initialize correlation vectors
cor.RT.k1.k2 <- rep(0,length(n.proj.seq))
cor.RT.maha <- rep(0,length(n.proj.seq))

  depth.RT.k2 <- depth.halfspace(x,x,num.directions = n.proj.k.star) # calculate RT with k star (= k2) random projections

for(i in 1:length(n.proj.seq)){
#  print(paste("n.proj.seq[i] =",n.proj.seq[i]))
#  print("k1")
  if(n.proj.seq[i] == 1){depth.RT.k1 <- depth.halfspace(x,x,num.directions = 2)}else{ # R function "depth.halfspace" does not work for k = 1 -_-')
  depth.RT.k1 <- depth.halfspace(x,x,num.directions = n.proj.seq[i]) # calculate RT with k1 random projections
  }
  cor.RT.k1.k2[i] <- cor(depth.RT.k1,depth.RT.k2)  # k2 is now equal to k.star
#  cor.RT.maha[i] <- cor(depth.RT.k1, depth.x.M.true)
}

# print(proc.time() - ptm)


if(length(cor.RT.k1.k2) == length(cor.RT.maha) & length(cor.RT.maha) == length(n.proj.seq)){
  data.correlation[1,1:length(n.proj.seq),(l-1)*length(dim.mv.simulations)+m] <- n.proj.seq
  data.correlation[2,1:length(n.proj.seq),(l-1)*length(dim.mv.simulations)+m] <- cor.RT.k1.k2
  data.correlation[3,1:length(n.proj.seq),(l-1)*length(dim.mv.simulations)+m] <- cor.RT.maha
}

}#end (for l in )
}#end (for m in)


}# end if(1==1) #run whole simulation 
# system("say Your R Code just finished!")

data.correlation.MC[,,,z] <- data.correlation

}



data.correlation.MC.mean <- apply(data.correlation.MC,1:3,mean) # mean over all mc simulations
dimnames(data.correlation.MC.mean)[[1]] <- c("k","RT","MH")

data.correlation.MC.05.quantile <- apply(data.correlation.MC,1:3, function(x) quantile(x, 0.05)) # 5-th quantile over all mc simulations
dimnames(data.correlation.MC.05.quantile)[[1]] <- c("k","RT","MH")

data.correlation.MC.95.quantile <- apply(data.correlation.MC,1:3, function(x) quantile(x, 0.95)) # 95-th quantile over all mc simulations
dimnames(data.correlation.MC.95.quantile)[[1]] <- c("k","RT","MH")

rho.to.reach <- 0.99
data.correlation.MC.mean.over.rho.to.reach <- data.correlation.MC.mean[2,,1] > rho.to.reach # correlation for RT.k vs RT.k.star over rho.to.reach?
data.correlation.MC.mean.first.to.reach.index <- min(which(data.correlation.MC.mean.over.rho.to.reach)) # first index to reach rho.to.reach
data.correlation.MC.mean.first.to.reach.k <- data.correlation.MC.mean[1,data.correlation.MC.mean.first.to.reach.index,1]

data.correlation.MC.05.quantile.over.rho.to.reach <- data.correlation.MC.05.quantile[2,,1] > rho.to.reach # correlation for RT.k vs RT.k.star over rho.to.reach?
data.correlation.MC.05.quantile.first.to.reach.index <- min(which(data.correlation.MC.05.quantile.over.rho.to.reach)) # first index to reach rho.to.reach
data.correlation.MC.05.quantile.first.to.reach.k <- data.correlation.MC.05.quantile[1,data.correlation.MC.05.quantile.first.to.reach.index,1]

data.correlation.MC.95.quantile.over.rho.to.reach <- data.correlation.MC.95.quantile[2,,1] > rho.to.reach # correlation for RT.k vs RT.k.star over rho.to.reach?
data.correlation.MC.95.quantile.first.to.reach.index <- min(which(data.correlation.MC.95.quantile.over.rho.to.reach)) # first index to reach rho.to.reach
data.correlation.MC.95.quantile.first.to.reach.k <- data.correlation.MC.95.quantile[1,data.correlation.MC.95.quantile.first.to.reach.index,1]


# plot(data.correlation.MC.mean[1,,],data.correlation.MC.mean[2,,], type = "l", 
#      xlab = "k", ylim = range(data.correlation.MC.95.quantile[2,,],data.correlation.MC.05.quantile[2,,]),
#      ylab = "rho")
# lines(data.correlation.MC.05.quantile[1,,],data.correlation.MC.05.quantile[2,,], col = "red")
# lines(data.correlation.MC.95.quantile[1,,],data.correlation.MC.95.quantile[2,,], col = "red")
# abline(h = rho.to.reach, lty = 2)
# abline(v = data.correlation.MC.mean.first.to.reach.k, lty = 2)
# abline(v = data.correlation.MC.05.quantile.first.to.reach.k, lty = 2)
# abline(v = data.correlation.MC.95.quantile.first.to.reach.k, lty = 2)

# 
# plot(data.correlation.MC.mean[1,,],data.correlation.MC.mean[2,,], type = "l", 
#      xlab = "k", ylim = range(data.correlation.MC.mean[2,,],data.correlation.MC.05.quantile[2,,]),
#      ylab = "rho")
# lines(data.correlation.MC.05.quantile[1,,],data.correlation.MC.05.quantile[2,,], col = "red")
# abline(h = rho.to.reach, lty = 2)
# abline(v = data.correlation.MC.mean.first.to.reach.k, lty = 2)
# abline(v = data.correlation.MC.05.quantile.first.to.reach.k, lty = 2)

result.matrix.final[,,1,z2] <- data.correlation.MC.95.quantile[,,1]
result.matrix.final[,,2,z2] <- data.correlation.MC.mean[,,1]
result.matrix.final[,,3,z2] <- data.correlation.MC.05.quantile[,,1]

  
k.zero.final[,z2] <- c(data.correlation.MC.95.quantile.first.to.reach.k, data.correlation.MC.mean.first.to.reach.k, data.correlation.MC.05.quantile.first.to.reach.k)


print((proc.time() - ptm.z2)[3]/60)

} #end mega loop

print((proc.time() - ptm)[3]/60)

system("say Your R Code just finished!")

# uncomment to save

# saveRDS(result.matrix.final, file = "./MT_functions/R_Data_choice_of_k_RTD/result.matrix.final.part1.rds")
# saveRDS(k.zero.final, file = "./MT_functions/R_Data_choice_of_k_RTD/k.zero.final.part1.rds")

# uncomment to load

load.results <- readRDS("./MT_functions/R_Data_choice_of_k_RTD/result.matrix.final.part1.rds")
load.k.results <- readRDS("./MT_functions/R_Data_choice_of_k_RTD/k.zero.final.part1.rds")

# k.r <- round(load.k.results)
# xtable(k.r[,4:6], display = rep("d",4))

# PRINT PLOTS

print = T



if(print == T){pdf("r_tilde_k_k_star.pdf",width = 10, height = 4.7)}

par(mfrow = c(1,2), mar = c(4.2, 4.2, 2, 2) + 0.1)

i = 1
rho.to.reach = 0.99
k.intervall.i <- 1:40
y.range.plot <- c(0.85,1)
# y.range.plot <- range(data.correlation.MC.95.quantile[2,],data.correlation.MC.05.quantile[2,])

data.correlation.MC.mean <- load.results[,k.intervall.i,2,i]
data.correlation.MC.05.quantile <- load.results[,k.intervall.i,1,i] 
data.correlation.MC.95.quantile <- load.results[,k.intervall.i,3,i]
data.correlation.MC.mean.first.to.reach.k <- load.k.results[2,i]
data.correlation.MC.05.quantile.first.to.reach.k <- load.k.results[1,i]
data.correlation.MC.95.quantile.first.to.reach.k <- load.k.results[3,i]


plot(data.correlation.MC.mean[1,],data.correlation.MC.mean[2,], type = "l",
     xlab = "k", ylim = y.range.plot,
     ylab = expression(widetilde(r)["k,k"^"*"*",P"["n"]*",mean"]),
     main = "n = 50, K = 2", cex.main = 0.85, lwd = 2)
lines(data.correlation.MC.05.quantile[1,],data.correlation.MC.05.quantile[2,], col = "darkgrey")
lines(data.correlation.MC.95.quantile[1,],data.correlation.MC.95.quantile[2,], col = "darkgrey")
abline(h = rho.to.reach, lty = 2)
abline(h = 1, lty = 2)
abline(v = data.correlation.MC.mean.first.to.reach.k, lty = 2)
abline(v = data.correlation.MC.05.quantile.first.to.reach.k, lty = 2)
abline(v = data.correlation.MC.95.quantile.first.to.reach.k, lty = 2)
text(data.correlation.MC.mean.first.to.reach.k+0.10*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,mean"]))
text(data.correlation.MC.95.quantile.first.to.reach.k+0.09*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,0.95"]))
text(data.correlation.MC.05.quantile.first.to.reach.k+0.07*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,0.05"]))
text(0.90*diff(range(k.intervall.i)),rho.to.reach - 0.05*(diff(y.range.plot)), label = expression(rho["0"]*" = 0.99"))

i = 4
rho.to.reach = 0.99
k.intervall.i <- 1:800
y.range.plot <- c(0.95,1)
# y.range.plot <- range(data.correlation.MC.95.quantile[2,],data.correlation.MC.05.quantile[2,])

data.correlation.MC.mean <- load.results[,k.intervall.i,2,i]
data.correlation.MC.05.quantile <- load.results[,k.intervall.i,1,i] 
data.correlation.MC.95.quantile <- load.results[,k.intervall.i,3,i]
data.correlation.MC.mean.first.to.reach.k <- load.k.results[2,i]
data.correlation.MC.05.quantile.first.to.reach.k <- load.k.results[1,i]
data.correlation.MC.95.quantile.first.to.reach.k <- load.k.results[3,i]


plot(data.correlation.MC.mean[1,],data.correlation.MC.mean[2,], type = "l",
     xlab = "k", ylim = y.range.plot,
     ylab = expression(widetilde(r)["k,k"^"*"*",P"["n"]*",mean"]),
     main = "n = 50, K = 3", cex.main = 0.85, lwd = 2)
lines(data.correlation.MC.05.quantile[1,],data.correlation.MC.05.quantile[2,], col = "darkgrey")
lines(data.correlation.MC.95.quantile[1,],data.correlation.MC.95.quantile[2,], col = "darkgrey")
abline(h = rho.to.reach, lty = 2)
abline(h = 1, lty = 2)
abline(v = data.correlation.MC.mean.first.to.reach.k, lty = 2)
abline(v = data.correlation.MC.05.quantile.first.to.reach.k, lty = 2)
abline(v = data.correlation.MC.95.quantile.first.to.reach.k, lty = 2)
text(data.correlation.MC.mean.first.to.reach.k+0.10*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,mean"]))
text(data.correlation.MC.95.quantile.first.to.reach.k+0.09*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,0.95"]))
text(data.correlation.MC.05.quantile.first.to.reach.k+0.06*diff(range(k.intervall.i)), 
     y.range.plot[1]+0.15*(diff(y.range.plot)), 
     label = expression("k"["0,0.05"]))
text(0.92*diff(range(k.intervall.i)),rho.to.reach - 0.05*(diff(y.range.plot)), label = expression(rho["0"]*" = 0.99"))


if(print == T){dev.off()}


