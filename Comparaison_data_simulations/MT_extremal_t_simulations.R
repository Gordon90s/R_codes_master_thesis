# set appropriate directory

main.directory <- "/Users/gordonschucker/Dropbox/MASTER/Masterarbeit & Paper/R & LaTeX"
setwd(main.directory)

library(SpatialExtremes) # to simulate max stable processes

#--------------------------------------------------------

ptm.et <- proc.time()
# Monte Carlo simulations of Extremal-t data
# saveRDS(BR.array, "BR_array_1000_500_2.rds")

set.seed(1)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 250 # Number of values for each stochastic process
mc <- 1000 # Number of Monte Carlo simulations

x.star <- 1:n.values
extremal.t.array <- array(NA, dim = c(n.draws, n.values, mc))
extremal.t.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
extremal.t.data <- rmaxstab(n.draws, x.star, "twhitmat", DoF = 4, nugget = 0, range = 3,
                            smooth = 0.7)
extremal.t.array[,,i] <- extremal.t.data

}                  

saveRDS(extremal.t.array, "./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_250_1000.rds")

time.et <- proc.time() - ptm.et
time.et

# ========================================================================
# ========================================================================
# ========================================================================


set.seed(2)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 100 # Number of values for each stochastic process
mc <- 1000 # Number of Monte Carlo simulations

x.star <- 1:n.values
extremal.t.array <- array(NA, dim = c(n.draws, n.values, mc))
extremal.t.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
extremal.t.data <- rmaxstab(n.draws, x.star, "twhitmat", DoF = 4, nugget = 0, range = 3,
                            smooth = 0.7)
extremal.t.array[,,i] <- extremal.t.data

}                  

saveRDS(extremal.t.array, "./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_100_1000.rds")

time.et <- proc.time() - ptm.et
time.et


# ========================================================================
# ========================================================================
# ========================================================================


set.seed(3)

n.draws <- 1000 # Number of stochastic process drawn
n.values <- 50 # Number of values for each stochastic process
mc <- 1000 # Number of Monte Carlo simulations

x.star <- 1:n.values
extremal.t.array <- array(NA, dim = c(n.draws, n.values, mc))
extremal.t.data <- matrix(rep(0,n.values*n.draws),nrow=n.draws)                  
                  
for(i in 1:mc){
extremal.t.data <- rmaxstab(n.draws, x.star, "twhitmat", DoF = 4, nugget = 0, range = 3,
                            smooth = 0.7)
extremal.t.array[,,i] <- extremal.t.data

}                  

saveRDS(extremal.t.array, "./MT_functions/Comparaison_data_simulations/extremal_t_array.rds_1000_50_1000.rds")

time.et <- proc.time() - ptm.et
time.et


# for testing purposes
# n.draws <- 2# Number of stochastic process drawn
# n.values <- 5 # Number of values for each stochastic process
# mc <- 2 # Number of Monte Carlo simulatio