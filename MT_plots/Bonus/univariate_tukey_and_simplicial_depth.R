# plot F(1-F) vs density
#---------------------------


library(Rlab)

x.norm <- seq(-3.5,3.5, by=0.01)
p.norm <- pnorm(x.norm)
d.norm <- dnorm(x.norm)

p.1.minus.F <- rep(0, length(x.norm))
for (i in 1:length(x.norm)){
  p.1.minus.F[i] <- min(1-p.norm[i],p.norm[i])
}

simplicial.p.1.minus.F <- rep(0, length(x.norm))
for (i in 1:length(x.norm)){
  simplicial.p.1.minus.F[i] <- (1-p.norm[i])*p.norm[i]
}

ed.F <- rep(0, length(x.norm))
for (i in 1:length(x.norm)){
  ed.F[i] <- 1 - abs((1-p.norm[i])-p.norm[i])
}



plot(x.norm, p.norm, type = "l", ylim = c(0,1), ylab = "y")
lines(x.norm, d.norm, col = "red")
lines(x.norm, p.1.minus.F, col = "blue")
legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))"), col = c("black","red","blue"), lty = c(1,1,1))

pdf("univariate_tukey_and_simplicial_depth.pdf", width = 9, height = 9)
par(mar = c(5,4,4,2)+0.1)
plot(x.norm, p.norm, type = "l", ylim = c(0,1), ylab = "y", xlab = "x")
lines(x.norm, d.norm, col = "red")
lines(x.norm, p.1.minus.F, col = "blue")
lines(x.norm, simplicial.p.1.minus.F, col = "blue", lty = 2)
legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))","(1-F(x))*F(x)"), 
       col = c("black","red","blue","blue"), lty = c(1,1,1,2))
dev.off()


plot(x.norm, p.norm, type = "l", ylim = c(0,1), ylab = "y")
lines(x.norm, d.norm, col = "red")
lines(x.norm, p.1.minus.F, col = "blue")
lines(x.norm, simplicial.p.1.minus.F, col = "blue", lty = 2)
lines(x.norm, 2*simplicial.p.1.minus.F, col = "blue", lty = 3)
lines(x.norm, ed.F, col = "darkred", lty = 2)
legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))","(1-F(x)*F(x)", "2*(1-F(x)*F(x)", "1-|F(x) - (1-F(x)-)|"), 
       col = c("black","red","blue","blue","blue","darkred"), lty = c(1,1,1,2,3,2))




#------------------------------------------------------

par(mar = c(5, 4, 4, 2) + 0.1)

set.seed(1)
x <- rgamma(1000)
plot(density(x))

range.x <- seq(0,10, length = 512)
depth.S <- univariate.simplicial.depth(range.x,x)
depth.T <- univariate.tukey.depth(range.x,x)

plot(range.x,depth.S*2, type = "l", col = "red")
lines(range.x,depth.T)







