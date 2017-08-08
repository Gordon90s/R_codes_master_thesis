# Master Thesis Plots

# Random Tukey Depth

#install.packages("mgcv")
library(mgcv)


print = T

set.seed(1)
V = matrix(c(1,0,0,1), nc = 2)
a <- rmvn(100, mu = c(0,0), V)
a.point <- c(1.03302370,  1.480399899)
p1 <- c(0.8064921,  0.59124484)
p2 <- c(0.69473088, -0.71926977)
ratio <- p1[2]/p1[1]

a.p1 <- p1 %*% t(a)
a.p1.density.x <- density(a.p1)$x
a.p1.density.y <- density(a.p1)$y
#a.p1.density.y <- rep(0, length(a.p1.density.x))

seq <- seq(-10,10,by = 0.01)
seq.p1 <- seq*ratio

scale <- 15
a.p1.density.y.transf <- a.p1.density.y*scale

graph.dim <- 1.1

if(print == T){
pdf("fig3.pdf", width = 10, height = 4.7)
}
par(mfrow = c(1,2), mar = c(4.2, 4, 1, 2) + 0.1, cex = 0.95)

lim <- 6
cex = 1/2
des.pol <- 40


seq.point <- -seq/ratio + a.point[2]


seq.polygon1 <- seq.polygon2 <- 355:512


polygon1 <- t(rbind(a.p1.density.x,a.p1.density.x*ratio))
polygon2 <- t(rbind(a.p1.density.x - a.p1.density.y.transf*p1[2],a.p1.density.x*ratio + a.p1.density.y.transf*p1[1]))
pol1.seq <- polygon1[seq.polygon1,]
pol2.seq <- polygon2[seq.polygon2,]
rev.pol2.seq <- polygon2[rev(seq.polygon2),]



plot(a, xlab = expression("x"[1]), ylab = expression("x"[2]), asp = 1, cex = cex, xlim = c(-lim,lim), ylim = c(-lim,lim))
polygon(rbind(pol1.seq,rev.pol2.seq), col = NA, angle = 80, border = NA, lty = 1,density = des.pol)
lines((a.p1.density.x - a.p1.density.y.transf*p1[2]),a.p1.density.x*ratio + a.p1.density.y.transf*p1[1], col = "red", lwd = 2)
lines(seq,seq.p1)
lines(seq+a.point[1], seq.point, col = "red", lty = 2)
points(a.point[1],a.point[2] , col = "red", cex = cex*3, pch = 16)
arrows(3,5,a.point[1]+0.3,a.point[2]+0.3, length = 0.15)
text(3.15,5.3, label = expression("min(F"["v"[1]]*"(x),1-F"["v"[1]]*"(x))"), col = "black")
text(-4.7,-4.2, label = expression("v"[1]))

arrows(5.3,0,a.point[1]+0.24,a.point[2]-0.03, length = 0.15, col = "red")
text(5.55,-0.05,label = "x", col = "red")


p1 <- c(0.69473088, -0.71926977)
ratio <- p1[2]/p1[1]

a.p1 <- p1 %*% t(a)
a.p1.density.x <- density(a.p1)$x
a.p1.density.y <- density(a.p1)$y
#a.p1.density.y <- rep(0, length(a.p1.density.x))

seq <- seq(-10,10,by = 0.01)
seq.p1 <- seq*ratio

scale <- 15
a.p1.density.y.transf <- a.p1.density.y*scale

lim <- 6
cex = 1/2

seq.point <- -seq/ratio + a.point[2]



seq.polygon1 <- seq.polygon2 <- 1:253


polygon1 <- t(rbind(a.p1.density.x,a.p1.density.x*ratio))
polygon2 <- t(rbind(a.p1.density.x - a.p1.density.y.transf*p1[2],a.p1.density.x*ratio + a.p1.density.y.transf*p1[1]))
pol1.seq <- polygon1[seq.polygon1,]
pol2.seq <- polygon2[seq.polygon2,]
rev.pol2.seq <- polygon2[rev(seq.polygon2),]

plot(a, xlab = expression("x"[1]), ylab = expression("x"[2]), asp = 1, cex = cex, xlim = c(-lim,lim), ylim = c(-lim,lim))
polygon(rbind(pol1.seq,rev.pol2.seq), col = NA, angle = 10, border = NA, lty = 1,density = des.pol)
lines((a.p1.density.x - a.p1.density.y.transf*p1[2]),a.p1.density.x*ratio + a.p1.density.y.transf*p1[1], col = "red", lwd = 2)

lines(seq,seq.p1)
lines(seq+a.point[1], seq.point, col = "red", lty = 2)
points(a.point[1],a.point[2] , col = "red", cex = cex*3, pch = 16)

arrows(0,5,a.point[1]-1.2,a.point[2]+0.8, length = 0.15)
text(-0.07,5.4,label = expression("min(F"["v"[2]]*"(x),1-F"["v"[2]]*"(x))"))
text(4.0,-4.85, label = expression("v"[2]))

arrows(5.3,0,a.point[1]+0.29,a.point[2]-0.03, length = 0.15, col = "red")
text(5.55,-0.05,label = "x", col = "red")

if(print == T){
dev.off()
}


par(mfrow = c(1,1))
