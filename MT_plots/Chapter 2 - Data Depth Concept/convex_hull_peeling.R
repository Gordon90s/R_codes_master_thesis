# Master Thesis Plots

# Random Tukey Depth

print = T

set.seed(20)

num.rv <- 25

a <- rmvn(num.rv, mu = c(0,0), V = matrix(c(1,0,0,1), nc = 2))

depth.a <- rep(0,num.rv)

for(i in 1:num.rv){
  depth.a[i] <- depth(a[i,],a)
}

depth.values <- c(0.05,0.1,0.20)
index.a1 <- depth.a > depth.values[1]
a1 <- a[index.a1,]
polygon.a <- a[chull(a),]
polygon.a1 <- a1[chull(a1),]

index.a2 <- depth.a > depth.values[2]
a2 <- a[index.a2,]
polygon.a <- a[chull(a),]
polygon.a2 <- a2[chull(a2),]


index.a3 <- depth.a > depth.values[3]
a3 <- a[index.a3,]
polygon.a <- a[chull(a),]
polygon.a3 <- a3[chull(a3),]


cex.fig <- 1.65

if(print == T){
pdf("fig1.pdf", width = 9, height = 9)
}
op <- par(mfrow = c(2,2),
          oma = c(0,0,0,0) + 0.5,
          mar = c(0,0,0,0) + 0.0)

plot(a, xlim = c(-3.3,2.2), ylim = c(-2.2,1.75), xaxt = "n", yaxt = "n", sub = "test")
polygon(polygon.a)
text(-0.4,-2.15,label = "polygon-0", cex = cex.fig)


plot(a, ann = F, xlim = c(-3.3,2.2), ylim = c(-2.2,1.75), xaxt = "n", yaxt = "n")
polygon(polygon.a)
polygon(polygon.a1)
text(-0.4,-2.15,label = "polygon-0 and 1", cex = cex.fig)

plot(a, ann = F, xlim = c(-3.3,2.2), ylim = c(-2.2,1.75), xaxt = "n", yaxt = "n")
polygon(polygon.a)
polygon(polygon.a1)
polygon(polygon.a2)
text(-0.4,-2.15,label = "polygon-0,1 and 2", cex = cex.fig)


plot(a, ann = F, xlim = c(-3.3,2.2), ylim = c(-2.2,1.75), xaxt = "n", yaxt = "n")
polygon(polygon.a)
polygon(polygon.a1)
polygon(polygon.a2)
polygon(polygon.a3)
text(-0.4,-2.15,label = "polygon-0,1,2 and 3", cex = cex.fig)

par(op)



if(print == T){
dev.off()
}


par(mfrow = c(1,1))
