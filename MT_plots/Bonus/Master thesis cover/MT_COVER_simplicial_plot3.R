# simplicial plot using jitter technique

library(mvtnorm)
library(scales)
library(fMultivar)
library(base)
library(ddalpha)


print = T

n.x <- 100

set.seed(4)

x.norm <- rmvnorm(n.x, mean = c(0,0)) 
#+ rcauchy2d(n.x, rho = 0)/50
jitter.factor <- 0
cex <- 1


depth.x <- rep(0,n.x)
depth.x <- depth.simplicial(x.norm,x.norm)
index <- 0
pair.i.j <- c(NA,NA)
depth.weight <- c(NA,NA)


grey.begin <- 0.0
grey.end <- 1
alpha.value <- 0.5
square.index <- 1.3


for(i in 1:n.x){
  for(j in 1:n.x){
    if(i < j){
      index <- index + 1
      if(index == 1){
        pair.i.j <- c(i,j)
        depth.weight <- (depth.x[i]+depth.x[j])^square.index
        }else{
        pair.i.j <- cbind(pair.i.j,c(i,j))
        depth.weight <- cbind(depth.weight,(depth.x[i]+depth.x[j])^square.index)
      }
    }
  }
}


depth.weight.grey <- abs((depth.weight-min(depth.weight))/(max(depth.weight-min(depth.weight)))-1)*(grey.end-grey.begin)+grey.begin
depth.weight.grey <- normalize(depth.weight)
order.depth.weight.grey <- order(depth.weight.grey, decreasing = T)


alpha.value.vector <- seq(1,1,length = index)

# col <- viridis_pal(alpha = alpha.value, option = "D", begin = grey.begin, end = grey.end)(index)
col.depth.weight <- rep(0,index)
for(i in 1:index){
  col.depth.weight[i] <- viridis_pal(alpha = alpha.value.vector[i], option = "D", begin = depth.weight.grey[i], end = depth.weight.grey[i], direction = 1)(1)
  if((i %% 10) == 0){print(i)}
}



line <- array(NA, dim = c(index,2,2))

#par(mar = c(6,5,5,2)+0.1)

if(print == T){pdf("MT_COVER_simplicial_depth.pdf", width = 3, height = 3)}
par(mar = c(0,0,0,0)+0.1)
plot(x.norm, pch = ".", ylab="", xlab="", main="", xaxt="n", yaxt="n") 
for(i in 1:index){
  line[i,,] <- cbind(x.norm[pair.i.j[1,i],],x.norm[pair.i.j[2,i],])
}

for(i in index:1){
  polygon(line[order.depth.weight.grey[i],1,],line[order.depth.weight.grey[i],2,], border = col.depth.weight[order.depth.weight.grey[i]])
}
points(x.norm, pch = ".")

if(print == T){dev.off()}
