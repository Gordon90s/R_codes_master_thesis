# set appropriate directory

main.directory <- "/Users/gordonschucker/Dropbox/MASTER/Masterarbeit & Paper/R & LaTeX"
setwd(main.directory)

library(SpatialExtremes) # to simulate max stable processes

#--------------------------------------------------------

ptm.et <- proc.time()

set.seed(1)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 25 # Number of values for each stochastic process
mc <- 1000 # Number of Monte Carlo simulations

x.star <- 1:n.values
brown.resnick.array <- array(NA, dim = c(n.draws, n.values, mc))
brown.resnick.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
brown.resnick.data <- rmaxstab(n.draws, x.star, "brown", range = 3, smooth = 0.7)
brown.resnick.array[,,i] <- brown.resnick.data

}                  

saveRDS(brown.resnick.array, "./MT_functions/Comparaison_data_simulations/BR_array_1000_25_1000.rds")

time.et <- proc.time() - ptm.et
time.et




ptm.et <- proc.time()

set.seed(2)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 50 # Number of values for each stochastic process
mc <- 1000 # Number of Monte Carlo simulations

x.star <- 1:n.values
brown.resnick.array <- array(NA, dim = c(n.draws, n.values, mc))
brown.resnick.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
brown.resnick.data <- rmaxstab(n.draws, x.star, "brown", range = 3, smooth = 0.7)
brown.resnick.array[,,i] <- brown.resnick.data

}                  

saveRDS(brown.resnick.array, "./MT_functions/Comparaison_data_simulations/BR_array_1000_50_1000.rds")

time.et <- proc.time() - ptm.et
time.et



ptm.et <- proc.time()

set.seed(3)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 100 # Number of values for each stochastic process
mc <- 10 # Number of Monte Carlo simulations

x.star <- 1:n.values
brown.resnick.array <- array(NA, dim = c(n.draws, n.values, mc))
brown.resnick.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
brown.resnick.data <- rmaxstab(n.draws, x.star, "brown", range = 3, smooth = 0.7)
brown.resnick.array[,,i] <- brown.resnick.data

}                  

saveRDS(brown.resnick.array, "./MT_functions/Comparaison_data_simulations/BR_array_1000_100_1000.rds")

time.et <- proc.time() - ptm.et
time.et

