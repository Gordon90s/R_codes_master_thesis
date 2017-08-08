# simulation of bivariate standard normal data colored by random Tukey depth

library(mvtnorm)
library(scales)
library(fMultivar)
library(base)
library(ddalpha)

print = T

n.x <- 2000

set.seed(4)

x.norm <- rmvnorm(n.x, mean = c(0,0)) 


depth.x.norm <- depth.halfspace(x.norm, x.norm, num.directions = 1000)
normalize.depth.x <- normalize(depth.x.norm)
order.depth.x <- order(depth.x.norm)

sort.x.norm <- x.norm[order.depth.x,]
normalized.depth.order <- normalize.depth.x[order.depth.x]


grey.begin <- 0.0
grey.end <- 1
# col <- viridis_pal(alpha = alpha.value.vector[i], option = "B", begin = grey.begin, end = grey.end)(n.x)
col <- viridis_pal(alpha = 1, option = "A", begin = grey.begin, end = grey.end)(n.x)[as.numeric(cut(normalized.depth.order, breaks = n.x))]


if(print == T){pdf("MT_COVER_random_tukey_vanilla.pdf", width = 3, height = 3)}

par(mfrow = c(1,1), mar = c(0,0,0,0)+0.1, cex = 1.6)

plot(x.norm, ylab="", xlab="", main="", xaxt="n", yaxt="n", col = "white")
for(i in 1:n.x){
  points(sort.x.norm[i,1],sort.x.norm[i,2], pch = "*", col = col[i])
}

if(print == T){dev.off()}
