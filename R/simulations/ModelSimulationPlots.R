library("MASS")
library("fields")
library("dplyr")
source("R/ModelSimulator.R")

pdf(file="Figures/Illustration.pdf",height=7.4,width=8)
par(mfrow=c(3,3)) 

set.seed(198345)

d <- 6 # number of parameters

theta_scenarios <- matrix(
  c("kappa" = 3,   "sigma" = 1,  "xi" = 0.05, "thL" = 10, "thU" = 10, "thw" = 0.25,  
    "kappa" = 0.3, "sigma" = 1,  "xi" = 0.05, "thL" = 10, "thU" = 0.5, "thw" = 0.25,  
    "kappa" = 3,   "sigma" = 1,  "xi" = 0.2,  "thL" = 4,  "thU" = 0.5, "thw" = 0.25),  
  nrow = d,  
  dimnames = list(c("kappa", "sigma", "xi", "thL", "thU", "thw"), NULL)
)

int_path <- "intermediates"
dir.create(int_path, recursive = TRUE, showWarnings = FALSE)
theta_scenarios %>% saveRDS(file = file.path(int_path, "theta_scenarios.rds"))

n <- 2000

### SIMULATION 1
theta <- theta_scenarios[,1]

sim <- rsim(n,theta)
Y <- sim$Y
Y1 <- Y[,1]
Y2 <- Y[,2]
R <- sim$R

ind.upp <- which(R>quantile(R,0.95))
ind.low <- which(R<quantile(R,0.05))

sim.big <- rsim(10^6,theta)
Y.big <- sim.big$Y
Y1.unif.big <- rank(Y.big[,1])/(10^6+1)
Y2.unif.big <- rank(Y.big[,2])/(10^6+1)
res <- kde2d(Y1.unif.big,Y2.unif.big,n=100)

par(mar=c(4,4,2,1),mgp=c(2.5,1,0))

plot(Y1,Y2,xlim=range(c(Y1,Y2)),ylim=range(c(Y1,Y2)),pch=20,xlab=expression(Y[1]),ylab=expression(Y[2]),asp=1,cex=0.5)
points(Y1[ind.upp],Y2[ind.upp],pch=20,col="red",cex=0.7)
points(Y1[ind.low],Y2[ind.low],pch=20,col="blue",cex=0.7)
abline(h=0,v=0,col="lightgrey")

plot(rank(Y1)/(n+1),rank(Y2)/(n+1),xlim=c(0,1),ylim=c(0,1),pch=20,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1,cex=0.5)
points((rank(Y1)/(n+1))[ind.upp],(rank(Y2)/(n+1))[ind.upp],pch=20,col="red",cex=0.7)
points((rank(Y1)/(n+1))[ind.low],(rank(Y2)/(n+1))[ind.low],pch=20,col="blue",cex=0.7)
abline(h=c(0,1),v=c(0,1),col="lightgrey")

par(mar=c(4,4,2,3),mgp=c(2.5,1,0))

image.plot(res$x,res$y,res$z,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1)

### SIMULATION 2
theta <- theta_scenarios[,2]

sim <- rsim(n,theta)
Y <- sim$Y
Y1 <- Y[,1]
Y2 <- Y[,2]
R <- sim$R

ind.upp <- which(R>quantile(R,0.95))
ind.low <- which(R<quantile(R,0.05))

sim.big <- rsim(10^6,theta)
Y.big <- sim.big$Y
Y1.unif.big <- rank(Y.big[,1])/(10^6+1)
Y2.unif.big <- rank(Y.big[,2])/(10^6+1)
res <- kde2d(Y1.unif.big,Y2.unif.big,n=100)

par(mar=c(4,4,2,1),mgp=c(2.5,1,0))

plot(Y1,Y2,xlim=range(c(Y1,Y2)),ylim=range(c(Y1,Y2)),pch=20,xlab=expression(Y[1]),ylab=expression(Y[2]),asp=1,cex=0.5)
points(Y1[ind.upp],Y2[ind.upp],pch=20,col="red",cex=0.7)
points(Y1[ind.low],Y2[ind.low],pch=20,col="blue",cex=0.7)
abline(h=0,v=0,col="lightgrey")

plot(rank(Y1)/(n+1),rank(Y2)/(n+1),xlim=c(0,1),ylim=c(0,1),pch=20,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1,cex=0.5)
points((rank(Y1)/(n+1))[ind.upp],(rank(Y2)/(n+1))[ind.upp],pch=20,col="red",cex=0.7)
points((rank(Y1)/(n+1))[ind.low],(rank(Y2)/(n+1))[ind.low],pch=20,col="blue",cex=0.7)
abline(h=c(0,1),v=c(0,1),col="lightgrey")

par(mar=c(4,4,2,3),mgp=c(2.5,1,0))

image.plot(res$x,res$y,res$z,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1)

### SIMULATION 3
theta <- theta_scenarios[,3]

sim <- rsim(n,theta)
Y <- sim$Y
Y1 <- Y[,1]
Y2 <- Y[,2]
R <- sim$R

ind.upp <- which(R>quantile(R,0.95))
ind.low <- which(R<quantile(R,0.05))

sim.big <- rsim(10^6,theta)
Y.big <- sim.big$Y
Y1.unif.big <- rank(Y.big[,1])/(10^6+1)
Y2.unif.big <- rank(Y.big[,2])/(10^6+1)
res <- kde2d(Y1.unif.big,Y2.unif.big,n=100)

par(mar=c(4,4,2,1),mgp=c(2.5,1,0))

plot(Y1,Y2,xlim=range(c(Y1,Y2)),ylim=range(c(Y1,Y2)),pch=20,xlab=expression(Y[1]),ylab=expression(Y[2]),asp=1,cex=0.5)
points(Y1[ind.upp],Y2[ind.upp],pch=20,col="red",cex=0.7)
points(Y1[ind.low],Y2[ind.low],pch=20,col="blue",cex=0.7)
abline(h=0,v=0,col="lightgrey")

plot(rank(Y1)/(n+1),rank(Y2)/(n+1),xlim=c(0,1),ylim=c(0,1),pch=20,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1,cex=0.5)
points((rank(Y1)/(n+1))[ind.upp],(rank(Y2)/(n+1))[ind.upp],pch=20,col="red",cex=0.7)
points((rank(Y1)/(n+1))[ind.low],(rank(Y2)/(n+1))[ind.low],pch=20,col="blue",cex=0.7)
abline(h=c(0,1),v=c(0,1),col="lightgrey")

par(mar=c(4,4,2,3),mgp=c(2.5,1,0))

image.plot(res$x,res$y,res$z,xlab=expression(Y[1]~"(on Unif(0,1) scale)"),ylab=expression(Y[2]~"(on Unif(0,1) scale)"),asp=1)

dev.off()

