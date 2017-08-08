# MT_MAIN_depth_function_for_analysis

# WARNING: run "Package_loader.R" first!



# this function loads all created functional depth function
# here we create specific "general depth function" that calculate multiple functional
# depths all at once for example


# general.int.depth(data)

# Input: 
#  1) data 
#     :: functional data set with which the depths are calculated
# Output: 
#     For
#        depth 1 : integrated depth with Tukey depth as marginal depth
#        depth 2 : integrated depth with Tukey depth as marginal depth and weight function w_3
#        depth 3 : integrated depth with simplicial depth as marginal depth
#        depth 4 : integrated depth with simplicial depth as marginal depth and weight function w_3
#        depth 5 : modified half region depth

#  1) depth in matrix form with nrow = 5 and ncol = nrow(data) (row 1 corresponds to depth 1, etc)
#  2) order in matrix form with nrow = 5 and ncol = nrow(data)
#  3) data.sorted in array form dim = c(nrow(data), ncol(data),5)

general.int.depth <- function(data){
  
number.of.depths <- 5  

depth.list <- vector("list", number.of.depths) 

depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

depth[1,] <- MDF(x = data, y = data, weight = "constant")$function.depth.x    
depth[3,] <- MDF(x = data, y = data)$function.depth.x  
depth[2,] <- MDF(x = data, y = data, weight = "constant", marginal.depth = "simplicial")$function.depth.x    
depth[4,] <- MDF(x = data, y = data, marginal.depth = "simplicial")$function.depth.x 
depth[5,] <- main.depth.function(data, graph = F, depth.type = 8)$depth


order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i])
}

return(depth.list)
}


# same as above, just different different depths. (used for creating various comparaisons)

#        depth 1 : integrated depth with Tukey depth as marginal depth
#        depth 2 : integrated depth with Tukey depth as marginal depth and weight function w_3
#        depth 3 : functional random projection depth with mean
#        depth 4 : functional random projection depth with min
#        depth 5 : extremal depth


general.depth.analysis <- function(data, n.proj = 1000, proj.type = "normal"){
  
number.of.depths <- 5  
depth.list <- vector("list", number.of.depths)   
  
depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

depth[1,] <- MDF(x = data, y = data, weight = "constant")$function.depth.x    
depth[2,] <- MDF(x = data, y = data)$function.depth.x  
depth.step1 <- MFRPD(data = data, more.than.one.depth = 2, n.proj = n.proj, proj.type = proj.type)
depth[3,] <- depth.step1$depth.mean    
depth[4,] <- depth.step1$depth.min 
depth[5,] <-  MDF(x = data, y = data, ED = T)$function.depth.x


order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i])
}

return(depth.list)
}












# same as above, just different different depths. (used for creating various comparaisons)


#        depth 1 : integrated depth with Tukey depth as marginal depth and weight function w_3
#        depth 2 : integrated depth with simplicial depth as marginal depth
#        depth 3 : integrated depth with simplicial depth as marginal depth and weight function w_3
#        depth 4 : extremal depth

general.for.outlier.depth <- function(data){
  
depth <- order <- matrix(NA, nrow = 4, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),4))  

depth[1,] <- MDF(x = data, y = data)$function.depth.x  
depth.step1 <- MFRPD(data = data, more.than.one.depth = 2, proj.type = "normal")
depth[2,] <- depth.step1$depth.mean    
depth[3,] <- depth.step1$depth.min     
depth[4,] <-  MDF(x = data, y = data, ED = T)$function.depth.x   


order <- t(apply(depth,1,order))

for(i in 1:4){
data.sorted[,,i] <- data[order[i,],]
}
return(list(depth = depth, order = order, data.sorted = data.sorted))

}

# same as above, but different depths

general.for.outlier.depth.EVT <- function(data, n.proj = 1000){
  
depth <- order <- matrix(NA, nrow = 3, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),3))  

depth.step1 <- MFRPD(data = data, more.than.one.depth = 4, EVT = T, n.proj = n.proj)
depth[1,] <- depth.step1$depth.mean
depth[2,] <- depth.step1$depth.min
depth[3,] <- depth.step1$depth.EVT.min


order <- t(apply(depth,1,order))

for(i in 1:3){
data.sorted[,,i] <- data[order[i,],]
}
return(list(depth = depth, order = order, data.sorted = data.sorted))

}



# to calculate the 10 most extremes functions for each depth

most.extremes.depth <- function(data, most.extremes = 10, col = iwanthue(nrow(data))){
  depth.data <- general.for.outlier.depth(data)
  
  
  col.outlier <- order.outlier <- depth.outlier <- matrix(NA, nrow = 4, ncol = most.extremes)  
  sorted.data.outlier <- array(NA, dim = c(most.extremes, ncol(data),4))
  
  order.outlier <- depth.data$order[,1:most.extremes]
  for(i in 1:4){
  depth.outlier[i,] <- depth.data$depth[i,order.outlier[i,]]
  col.outlier[i,] <- col[order.outlier[i,]]
  sorted.data.outlier[,,i] <- depth.data$data.sorted[1:most.extremes,,i]
  }
  
  return(list(depth.outlier = depth.outlier, order.outlier = order.outlier, 
              sorted.data.outlier = sorted.data.outlier, col.outlier = col.outlier, col = col))

}




# same as above, but with different depth functions

most.extremes.depth.EVT <- function(data, most.extremes = 10, col = iwanthue(nrow(data))){
  depth.data <- general.for.outlier.depth.EVT(data, n.proj = 1000)
  
  
  col.outlier <- order.outlier <- depth.outlier <- matrix(NA, nrow = 3, ncol = most.extremes)  
  sorted.data.outlier <- array(NA, dim = c(most.extremes, ncol(data),3))
  
  order.outlier <- depth.data$order[,1:most.extremes]
  for(i in 1:3){
  depth.outlier[i,] <- depth.data$depth[i,order.outlier[i,]]
  col.outlier[i,] <- col[order.outlier[i,]]
  sorted.data.outlier[,,i] <- depth.data$data.sorted[1:most.extremes,,i]
  }
  
  return(list(depth.outlier = depth.outlier, order.outlier = order.outlier, 
              sorted.data.outlier = sorted.data.outlier, col.outlier = col.outlier, col = col))

}



# BONUS functions, new depths

general.proj.depth <- function(data){
  
depth <- order <- matrix(NA, nrow = 2, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),2))  

depth.step1 <- MFRPD(data = data, more.than.one.depth = 2)

depth[1,] <- depth.step1$depth.mean    
depth[2,] <- depth.step1$depth.min  

order <- t(apply(depth,1,order))

for(i in 1:2){
data.sorted[,,i] <- data[order[i,],]
}
return(list(depth = depth, order = order, data.sorted = data.sorted))

}


# old function
# general.proj.EVT.depth2 <- function(data){
#   
# depth <- order <- matrix(NA, nrow = 4, ncol = nrow(data))
# data.sorted <- array(NA, dim = c(nrow(data), ncol(data),4))  
# 
# depth.step1 <- MFRPD(data = data, more.than.one.depth = 4, EVT = T)
# 
# depth[1,] <- depth.step1$depth.mean    
# depth[2,] <- depth.step1$depth.min  
# depth[3,] <- depth.step1$depth.EVT.mean    
# depth[4,] <- depth.step1$depth.EVT.min 
# 
# order <- t(apply(depth,1,order))
# 
# for(i in 1:4){
# data.sorted[,,i] <- data[order[i,],]
# }
# return(list(depth = depth, order = order, data.sorted = data.sorted))
# 
# }

# for analysis of projection depth with lists





# depth analysis for EVT for the FRPD


general.proj.EVT.depth <- function(data,n.proj = 1000){
number.of.depths <- 3 
depth.list <- vector("list", number.of.depths)   
  
depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

  
depth.step1 <- MFRPD(data = data, more.than.one.depth = 4, EVT = T, n.proj = n.proj)
depth[1,] <- depth.step1$depth.mean
depth[2,] <- depth.step1$depth.min
depth[3,] <- depth.step1$depth.EVT.min

   
# depth[6,] <- MDF(x = data, y = data, weight = "constant", marginal.depth = "simplicial")$function.depth.x    
# depth[7,] <- MDF(x = data, y = data, marginal.depth = "simplicial")$function.depth.x 


order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i], gamma.final = depth.step1$final.gamma)
}

return(depth.list)
}















#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------  NOT USED (?)------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++





# depth analysis for EVT
general.proj.EVT.RV.depth <- function(data,n.proj = 1000){
number.of.depths <- 6  
depth.list <- vector("list", number.of.depths)   
  
depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

  

depth.step1 <- MFRPD(data = data, more.than.one.depth = 6, EVT = T, n.proj = n.proj, proj.type = "abs")
  
depth[1,] <- depth.step1$depth.mean
depth[2,] <- depth.step1$depth.EVT.mean
depth[3,] <- depth.step1$depth.EVT.variable.gamma.mean
depth[4,] <- depth.step1$depth.min
depth[5,] <- depth.step1$depth.EVT.min
depth[6,] <- depth.step1$depth.EVT.variable.gamma.min


order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i])
}

return(depth.list)
}


# REAL ONE USED


# depth analysis for EVT only with "functional RV"
general.proj.EVT.FRV.depth <- function(data,n.proj = 1000){
number.of.depths <- 3  
depth.list <- vector("list", number.of.depths)   
  
depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

  

depth.step1 <- MFRPD(data = data, more.than.one.depth = 4, EVT = T, n.proj = n.proj, proj.type = "abs")
  
depth[1,] <- depth.step1$depth.mean
depth[2,] <- depth.step1$depth.min
depth[3,] <- depth.step1$depth.EVT.min



order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i])
}

return(depth.list)
}

























# NOT USED










# general.depth(data)

# Input: 
#  1) data 
#     :: functional data set with which the depths are calculated
# Output: 
# For
#        depth 1 : integrated depth with Tukey depth as marginal depth
#        depth 2 : integrated depth with Tukey depth as marginal depth and weight function w_3
#        depth 3 : functional random projection depth with Tukey depth as marginal depth, nproj = 1000 and mean function
#        depth 4 : functional random projection depth with Tukey depth as marginal depth, nproj = 1000 and min function
#        depth 5 : extremal depth

#  1) depth in matrix form with nrow = 5 and ncol = nrow(data) (row 1 corresponds to depth 1, etc)
#  2) order in matrix form with nrow = 5 and ncol = nrow(data)
#  3) data.sorted in array form dim = c(nrow(data), ncol(data),5)

general.depth <- function(data){
  
number.of.depths <- 5

depth.list <- vector("list", number.of.depths) 
    
depth <- order <- matrix(NA, nrow = number.of.depths, ncol = nrow(data))
data.sorted <- array(NA, dim = c(nrow(data), ncol(data),number.of.depths))  

depth[1,] <- MDF(x = data, y = data, weight = "constant")$function.depth.x    
depth[2,] <- MDF(x = data, y = data)$function.depth.x  
depth.step1 <- MFRPD(data = data, more.than.one.depth = 2)
depth[3,] <- depth.step1$depth.mean    
depth[4,] <- depth.step1$depth.min 
depth[5,] <-  MDF(x = data, y = data, ED = T)$function.depth.x   

order <- t(apply(depth,1,order))

for(i in 1:number.of.depths){
data.sorted[,,i] <- data[order[i,],]
}

# save as list
for(i in 1:number.of.depths){
  depth.list[[i]] <- list(depth = depth[i,], order = order[i,], data.sorted = data.sorted[,,i])
}

return(depth.list)

}






