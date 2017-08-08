# Choice of k FRPD

set.seed(2)

reload.packages = F

if(reload.packages == T){

#install.packages("fda.usc") # package for functional depths
library(fda.usc)
library(clusterGeneration)
#install.packages("xtable")
library(xtable) # to create LaTeX tables
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
}


  calcSigma <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
  for(i in 1:nrow(Sigma)){
  for (j in 1:ncol(Sigma)) Sigma[i,j]<-exp(-1/2*(abs(X1[i]-X2[j])/l)^2)
  }
  return(Sigma)
  }

# GOOD NEWS: for the functional random projeciton depth, the work is reproducible as the function depth.RT
  
n.draws.all <- c(50)
n.values.all <- c(50,100)
  
# n.draws.all <- c(100,250,1000)
# n.values.all <- c(50,100,250)

# n.proj.k.star <- 10
n.proj.k.star <- 10000



all.simulations <-  matrix(NA,length(n.draws.all)*length(n.values.all),2)
for(i in 1:length(n.draws.all)){
  for(j in 1:length(n.values.all)){
    all.simulations[j+(i-1)*length(n.values.all),] <- c(n.draws.all[i],n.values.all[j])
    }
  }



# n.proj <- 10
n.proj <- 3000


proj.seq <- seq(250,n.proj, length = n.proj/250)
n.proj.seq <- length(proj.seq)


result.matrix.final <- array(NA, dim = c(2,n.proj.seq,3,nrow(all.simulations)))
k.zero.final <- matrix(NA, 3, nrow(all.simulations))

ptm <- proc.time() # measure time of the code

for(z2 in 1:nrow(all.simulations)){
print(z2)
ptm.z2 <- proc.time()
# smaller test

# number of MC simulations
mc <- 1000


# simulations
n.values <- all.simulations[z2,2] # slightly change in the algorithm, keep length of 1. (only for naming though, so no big deal)
n.draws <- all.simulations[z2,1]

amount.of.simulations <- length(n.values)*length(n.draws)
data.correlation.MC <- array(NA, dim = c(2,n.proj.seq,amount.of.simulations,mc))
dimnames(data.correlation.MC)[[1]] <- c("k","FRPD")

for(z in 1:mc){
# print(z)
if(1==1){

data.correlation <- array(NA, dim = c(2,length(proj.seq),amount.of.simulations)) # to save the seq.sequences and the correlation values
dimnames(data.correlation)[[1]] <- c("k","cor FRPT.k.k.star")
    
  
# simulation Gaussian Process  
  
  
# The standard deviation of the noise
x.star <- seq(-5,5, len=n.values)
nval <- 5
f <- data.frame(x=seq(-5,5,l=nval), y=rnorm(nval,0,10))
sigma.n <- 1
# Recalculate the mean and covariance functions
k.xx <- calcSigma(f$x,f$x)
k.xxs <- calcSigma(f$x, x.star)
k.xsx <- calcSigma(x.star,f$x)
k.xsxs <- calcSigma(x.star,x.star)
f.bar.star <- k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs-k.xsx%*%solve(k.xx+sigma.n^2*diag(1,ncol(k.xx)))%*%k.xxs
values <- matrix (rep(0,length(x.star)*n.draws), ncol=n.draws)
for (i in 1:n.draws)  values[,i] <- mvrnorm(1,f.bar.star,cov.f.star)
smooth.gp.data <- t(values) # save data, need to transpose for having good matrix

# rainbow.plot(smooth.gp.data, x.star, col = "black")
# rainbow.plot(smooth.gp.data, x.star)

data <- smooth.gp.data


# initialize correlation vectors
cor.FRPD.k.k.star <- rep(0,length(proj.seq))

depth.FRPD.k.star <- depth.RT(data, nproj = n.proj.k.star)$dep # calculate FRTD with k.star random projections
  
for(i1 in 1:length(proj.seq)){
#  print(paste("n.proj.seq[i] =",n.proj.seq[i]))
#  print("k1")
  depth.FRPD.k <- depth.RT(data, nproj = proj.seq[i1])$dep # calculate FRTD with k random projections
  cor.FRPD.k.k.star[i1] <- cor(depth.FRPD.k,depth.FRPD.k.star)
}

data.correlation[1,,1] <- proj.seq
data.correlation[2,,1] <- cor.FRPD.k.k.star

# print(proc.time() - ptm)


}# end if(1==1) #run whole simulation 
# system("say Your R Code just finished!")

data.correlation.MC[,,,z] <- data.correlation

}



data.correlation.MC.mean <- apply(data.correlation.MC,1:3,mean) # mean over all mc simulations
dimnames(data.correlation.MC.mean)[[1]] <- c("k","FRPD")

data.correlation.MC.05.quantile <- apply(data.correlation.MC,1:3, function(x) quantile(x, 0.05)) # 5-th quantile over all mc simulations
dimnames(data.correlation.MC.05.quantile)[[1]] <- c("k","FRPD")

data.correlation.MC.95.quantile <- apply(data.correlation.MC,1:3, function(x) quantile(x, 0.95)) # 95-th quantile over all mc simulations
dimnames(data.correlation.MC.95.quantile)[[1]] <- c("k","FRPD")

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

# saveRDS(result.matrix.final, file = "./MT_functions/R_Data_choice_of_k_FRPD/result.matrix.final.FRPD.part2.rds")
# saveRDS(k.zero.final, file = "./MT_functions/R_Data_choice_of_k_FRPD/k.zero.final.FRPD.part2.rds")

# uncomment to load 

# readRDS("./MT_functions/R_Data_choice_of_k_FRPD/result.matrix.final.FRPD.part2.rds")
# readRDS("./MT_functions/R_Data_choice_of_k_FRPD/k.zero.final.FRPD.part2.rds")

k.zero.final
