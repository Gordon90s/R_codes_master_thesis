# Plot examples of scalars v_i, x_j


print = F

# load presimulated data
BR.data <- readRDS("./MT_functions/Comparaison_data_simulations/BR_array_1000_100_100.rds")
data <- BR.data[1:1000,,1]

n.draws <- dim(data)[1]
n.values <- dim(data)[2]

# =============================================================
# ------- MULTIVARIATE FUNCTIONAL PROJECTION DEPTH ------------
# =============================================================


# generate k projections

k <- 100
set.seed(10)

proj.functions <- rproc2fdata(k, t = 1:n.values)
proj <- proj.functions$data # regular Gaussian processes (as in functional random Tukey depth)

scalars <- array(0, dim = c(n.draws, k)) # calculate projections
for(i in 1:n.draws){
  for(j in 1:k){
    scalars[i,j] <- mean(data[i,]*proj[j,])
  }                      
}

# apply Tukey depth to projections
depth.k <- order.depth.k <- matrix(0, nrow = n.draws, ncol = k)
for(j in 1:k){
  depth.k[,j] <- univariate.tukey.depth(scalars[,j],scalars[,j]) # calculated tukey depth tukey depth (?)
  order.depth.k[,j] <- order(depth.k[,j])
}


most.extremes <- 5 # number of extrems to plot differently


# plot points, extremes are colored differently for each different k using depth calculated for each k

#--------------------------
# plot with density

k.1 <- 9 # number of plots

# par(mfrow = c(2,2), mar = c(5,4,4,3)+0.1)
par(mfrow = c(3,3), mar = c(4.2, 4, 1, 2))
points.data <- rep(0,n.draws)
for(j in 1:k.1){
  points.data <- scalars[order.depth.k[,j],j] # arrange data by depth for each different k
  plot(density(points.data),main = paste0("k = ",j), xlab = "y", ylab = "f(y)")
  points(points.data,rep(0,n.draws), pch = "*", cex = 1.2)
    for(l in 1:most.extremes){
    points(points.data[l],0, col = "red", cex = 2.5, pch = "*")
  }
}

# NEW PLOT WITH ONLY EXTREMES

k.1 <- c(1,3,4,9) # number of plots

# par(mfrow = c(2,2), mar = c(5,4,4,3)+0.1)

if(print == T){pdf("brown_resnick_selected_projections.pdf", height = 10, width = 10)}
par(mfrow = c(2,2), mar = c(4.2, 4, 1, 2))
points.data <- rep(0,n.draws)
for(j in k.1){
  points.data <- scalars[order.depth.k[,j],j] # arrange data by depth for each different k
  plot(density(points.data),main = paste0("i = ",j), xlab = expression("y"["i,n"]), ylab = expression("f(y"["i,n"]*")"))
  points(points.data,rep(0,n.draws), pch = "*", cex = 1.2)
    for(l in 1:most.extremes){
    points(points.data[l],0, col = "red", cex = 2.5, pch = "*")
  }
}
if(print == T){dev.off()}


