# simplicial plot using jitter technique

library(mvtnorm)
library(scales)
library(fMultivar)
library(base)
library(ddalpha)


print = F

n.x <- 100

set.seed(4)

x.norm <- rmvnorm(n.x, mean = c(0,0)) + rcauchy2d(n.x, rho = 0)/50
jitter.factor <- 0
cex <- 1


depth.x <- rep(0,n.x)
depth.x <- depth.simplicial(x.norm,x.norm)
index <- 0
pair.i.j <- c(NA,NA)
depth.weight <- c(NA,NA)


grey.begin <- 0.0
grey.end <- 0.95
alpha.value <- 0.5
square.index <- 1.75

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
order.depth.weight.grey <- order(depth.weight.grey, decreasing = T)



line <- array(NA, dim = c(index,2,2))

#par(mar = c(6,5,5,2)+0.1)
par(mar = c(5,4,4,1)+0.1)
if(print == T){pdf("simplicial_depth_2.pdf", width = 9, height = 9)}
plot(x.norm, pch = ".", xlab = expression("x"[1]), ylab = expression("x"[2]), cex.axis = cex, cex.lab = cex) 
for(i in 1:index){
  line[i,,] <- cbind(x.norm[pair.i.j[1,i],],x.norm[pair.i.j[2,i],])
}

for(i in 1:index){
  polygon(line[order.depth.weight.grey[i],1,],line[order.depth.weight.grey[i],2,], border = grey(depth.weight.grey[order.depth.weight.grey[i]], alpha = alpha.value))
}
points(x.norm, pch = ".")

if(print == T){dev.off()}
