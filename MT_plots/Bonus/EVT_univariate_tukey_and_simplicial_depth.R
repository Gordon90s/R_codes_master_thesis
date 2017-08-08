# plot F(1-F) vs density
#---------------------------

# SAME FOR FRECHET

x.norm <- seq(-1,25, by=0.01)
p.norm <- pgev(x.norm, shape = 1)
d.norm <- dgev(x.norm, shape = 1)

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



# plot(x.norm, p.norm, type = "l", ylim = c(0,1), ylab = "y")
# lines(x.norm, d.norm, col = "red")
# lines(x.norm, p.1.minus.F, col = "blue")
# legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))"), col = c("black","red","blue"), lty = c(1,1,1))

pdf("EVT_univariate_tukey_and_simplicial_depth_two_sided.pdf", width = 9, height = 9)
par(mar = c(5,4,4,2)+0.1)
plot(c(rev(-(x.norm +1)),x.norm +1), c(rev(-p.norm/2+0.5),p.norm/2+0.5), type = "l", ylim = c(0,1), ylab = "y", xlab = "x", xlim = c(-25,25))
lines(c(rev(-(x.norm +1)),x.norm +1), c(rev(d.norm/2),d.norm/2), col = "red")
lines(c(rev(-(x.norm +1)),x.norm +1), c(rev(p.1.minus.F),p.1.minus.F), col = "blue")
lines(c(rev(-(x.norm +1)),x.norm +1), c(rev(simplicial.p.1.minus.F),simplicial.p.1.minus.F), col = "blue", lty = 2)
legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))","(1-F(x))*F(x)"), 
       col = c("black","red","blue","blue"), lty = c(1,1,1,2))
dev.off()

pdf("EVT_univariate_tukey_and_simplicial_depth_one_sided.pdf", width = 9, height = 9)
par(mar = c(5,4,4,2)+0.1)
plot(x.norm +1, p.norm, type = "l", ylim = c(0,1), ylab = "y", xlab = "x", xlim = c(0,25))
lines(x.norm +1, d.norm, col = "red")
lines(x.norm +1, p.1.minus.F, col = "blue")
lines(x.norm +1, simplicial.p.1.minus.F, col = "blue", lty = 2)
legend(x = 13.0, y = 0.6, c("F(x)","f(x)","min(1-F(x),F(x))","(1-F(x))*F(x)"), 
       col = c("black","red","blue","blue"), lty = c(1,1,1,2))
dev.off()

# plot(x.norm, p.norm, type = "l", ylim = c(0,1), ylab = "y")
# lines(x.norm, d.norm, col = "red")
# lines(x.norm, p.1.minus.F, col = "blue")
# lines(x.norm, simplicial.p.1.minus.F, col = "blue", lty = 2)
# lines(x.norm, 2*simplicial.p.1.minus.F, col = "blue", lty = 3)
# lines(x.norm, ed.F, col = "darkred", lty = 2)
# legend("topleft", c("F(x)","f(x)","min(1-F(x),F(x))","(1-F(x)*F(x)", "2*(1-F(x)*F(x)", "1-|F(x) - (1-F(x)-)|"), 
#        col = c("black","red","blue","blue","blue","darkred"), lty = c(1,1,1,2,3,2))


test <- rgev(1000, shape = 1) #*rnorm(1000)
test <- (rgev(1000, shape = 1)+1) #*sample(c(-1,1), replace = T, size = 1000)
test.less <- sort(test)[-c(970:1000)]

plot(density(test.less), xlim = c(-10,10))
lines(density(-test.less), col = "red")

plot(density(test.less+1), xlim = c(-10,10))
lines(density(-test.less-1), col = "red")

plot(density(test.less-0.75), xlim = c(-10,10))
lines(density(-test.less+0.75), col = "red")
