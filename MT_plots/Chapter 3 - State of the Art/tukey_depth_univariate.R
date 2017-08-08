# Univariate tukey depth

print <- T


resolution <- 100
seq <- seq(-3.5,3.5,by= 1/resolution)
dnorm.seq <- dnorm(seq)
abline.v <- 0.45
des.pol <- 12
count.to.abline.v <- (abline.v - min(seq))*resolution



seq.count.to.abline.v <- 1:count.to.abline.v
seq.rest <- (count.to.abline.v+1):length(seq)

polygon1 <- t(rbind(seq[seq.count.to.abline.v],dnorm.seq[seq.count.to.abline.v]))
polygon2 <- t(rbind(seq[seq.rest],dnorm.seq[seq.rest]))
pol1.seq <- polygon1[seq.polygon1,]
pol2.seq <- polygon2[seq.polygon2,]
rev.pol1.seq <- polygon1[rev(seq.polygon1),]
rev.pol2.seq <- polygon2[rev(seq.polygon2),]

if(print == T){
pdf("tukey_depth_univariate.pdf", width = 9, height = 7.9)
}
par(mfrow = c(1,1), mar = c(4.2, 4, 1, 2) + 0.1, cex = 1.6)

plot(seq, dnorm.seq+0.003, type = "l", xlab = "x", ylab = "f(x)", ylim = c(0,0.55), lwd = 3, col = "red")
abline(v = abline.v, lty = 2, col = "red", lwd = 3)
polygon(x = c(seq[seq.count.to.abline.v], rev(seq[seq.count.to.abline.v])), y = c(rep(0,length(seq.count.to.abline.v)),rev(dnorm.seq[seq.count.to.abline.v])), col = NA, angle = 40, border = NA, lty = 1, density = des.pol, lwd = 2)
polygon(x = c(seq[seq.rest], rev(seq[seq.rest])), y = c(rep(0,length(seq.rest)),rev(dnorm.seq[seq.rest])), col = NA, angle = 130, border = NA, lty = 1, density = des.pol, lwd = 2)
axis(1,0.45,label = expression("x"[1]), col = "red", col.axis = "red")

#text(0.75,0.45, label = expression("x"[1]))


arrows(-1.84,0.40,-0.5,0.15, length = 0.15)
text(-1.990,0.43,label = expression("F(x"[1]*")"))


arrows(1.74,0.27,1,0.12, length = 0.15)
text(2.3,0.3,label = expression("1-F(x"[1]*") = D"["T"[1]]*"(x"[1]*", P)"))


if(print == T){
dev.off()
}
