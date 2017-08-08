# k RT vs k MH

reload.packages = F
save.data = F

if(reload.packages == T){

#install.packages("ddalpha")
library(ddalpha)
#install.packages("clusterGeneration")
library(clusterGeneration)
#install.packages("ddalpha")
library(xtable) # to create LaTeX tables
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
}


# WARNING: work not (exactly) reproducible as seed function within "depth.halfspace" is not working. (as of 12.04.2017)


dim.mv.simulations.all <- c(4)
n.x.simulations.all <- c(100)

# n.proj.k.star <- 10
n.proj.k.star <- 1000000



all.simulations <-  matrix(NA,length(dim.mv.simulations.all)*length(n.x.simulations.all),2)
for(i in 1:length(dim.mv.simulations.all)){
  for(j in 1:length(n.x.simulations.all)){
    all.simulations[j+(i-1)*length(n.x.simulations.all),] <- c(dim.mv.simulations.all[i],n.x.simulations.all[j])
    }
  }



# n.proj <- 10
n.proj <- 2100


n.proj.seq <- 1:n.proj

result.matrix.final <- array(NA, dim = c(3,n.proj,3,nrow(all.simulations)))
k.zero.final <- matrix(NA, 3, nrow(all.simulations))

ptm <- proc.time() # measure time of the code

z2 <- 1

ptm.z2 <- proc.time()
# smaller test

# number of MC simulations
mc <- 1


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
  cor.RT.maha[i] <- cor(depth.RT.k1, depth.x.M.true) 
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



}



print((proc.time() - ptm.z2)[3]/60)


print((proc.time() - ptm)[3]/60)

data.correlation <- adrop(data.correlation[,,1,drop=FALSE],drop=3)
if(save.data == T){saveRDS(data.correlation, file = "./MT_functions/R_Data_choice_of_k_RTD/data.correlation.MH4.rds")}

# PRINT PLOTS





dim.mv.simulations.all <- c(8)
n.x.simulations.all <- c(100)

# n.proj.k.star <- 10
n.proj.k.star <- 1000000



all.simulations <-  matrix(NA,length(dim.mv.simulations.all)*length(n.x.simulations.all),2)
for(i in 1:length(dim.mv.simulations.all)){
  for(j in 1:length(n.x.simulations.all)){
    all.simulations[j+(i-1)*length(n.x.simulations.all),] <- c(dim.mv.simulations.all[i],n.x.simulations.all[j])
    }
  }



# n.proj <- 10
n.proj <- 10500


n.proj.seq <- 1:n.proj

result.matrix.final <- array(NA, dim = c(3,n.proj,3,nrow(all.simulations)))
k.zero.final <- matrix(NA, 3, nrow(all.simulations))

ptm <- proc.time() # measure time of the code

z2 <- 1

ptm.z2 <- proc.time()
# smaller test

# number of MC simulations
mc <- 1


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
  cor.RT.maha[i] <- cor(depth.RT.k1, depth.x.M.true) 
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



}



print((proc.time() - ptm.z2)[3]/60)


print((proc.time() - ptm)[3]/60)

system("say Your R Code just finished!")


data.correlation <- adrop(data.correlation[,,1,drop=FALSE],drop=3)
if(save.data == T){saveRDS(data.correlation, file = "./MT_functions/R_Data_choice_of_k_RTD/data.correlation.MH8.rds")}


# PRINT PLOTS

print = T

data.correlation <- readRDS(file = "./MT_functions/R_Data_choice_of_k_RTD/data.correlation.MH4.rds")
n.proj <- ncol(data.correlation)
seq.plot <- seq(5,n.proj,length = n.proj/5)

if(print == T){pdf("r_tilde_k_MH.pdf",width = 10, height = 5)}
par(mfrow = c(1,2), mar = c(4.2, 4.2, 2, 2) + 0.1)

i = 1
rho.to.reach = 0.99
y.range.plot <- c(0.94,1)
# y.range.plot <- range(data.correlation[2,],data.correlation[3,])

plot(data.correlation[1,seq.plot],data.correlation[2,seq.plot], type = "l",
     xlab = "k", ylim = y.range.plot,
     ylab = expression(widetilde(r)["k,k"^"*"*",P"["n"]]*", "*"r"["k,MH,P"["n"]]),
     main = "n = 100, K = 4", cex.main = 0.85)
lines(data.correlation[1,seq.plot],data.correlation[3,seq.plot], col = "red")
abline(h = 1, lty = 2)
text(1300, 0.9805, label = expression(widetilde(r)["k,k"^"*"*",P"["n"]]))
arrows(1300, 0.982, 1200, 0.992, length = 0.13)
text(1250, 0.9505, label = expression("r"["k,MH,P"["n"]]), col = "red")
arrows(1250, 0.952, 1150, 0.962, length = 0.13, col = "red")



data.correlation <- readRDS(file = "./MT_functions/R_Data_choice_of_k_RTD/data.correlation.MH8.rds")
n.proj <- ncol(data.correlation)
seq.plot <- seq(20,n.proj,length = n.proj/20)

i = 1
rho.to.reach = 0.99
y.range.plot <- c(0.8,1)
# y.range.plot <- range(data.correlation[2,],data.correlation[3,])

plot(data.correlation[1,seq.plot],data.correlation[2,seq.plot], type = "l",
     xlab = "k", ylim = y.range.plot,
     ylab = expression(widetilde(r)["k,k"^"*"*",P"["n"]]*", "*"r"["k,MH,P"["n"]]),
     main = "n = 100, K = 8", cex.main = 0.85)
lines(data.correlation[1,seq.plot],data.correlation[3,seq.plot], col = "red")
abline(h = 1, lty = 2)
text(1700, 0.9805, label = expression(widetilde(r)["k,k"^"*"*",P"["n"]]))
arrows(1850, 0.975, 2900, 0.97, length = 0.13)
text(6250, 0.84, label = expression("r"["k,MH,P"["n"]]), col = "red")
arrows(6250, 0.845, 6000, 0.875, length = 0.13, col = "red")

if(print == T){dev.off()}


