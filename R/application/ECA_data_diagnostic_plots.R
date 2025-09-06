##################################################
##  Pairs plots btw the real and simulated data ##
##################################################
parameter_labels <- c(
  "kappa" = expression(kappa), 
  "sigma" = expression(sigma), 
  "xi" = expression(xi), 
  "thL" = expression(theta[L]), 
  "thU" = expression(theta[U]), 
  "thw" = expression(theta[omega])
) 
set.seed(123)

n11 <- ncol(mat_11_15_scaled) #2138
n22 <- ncol(mat_11_23_scaled) #2166
n33 <- ncol(mat_15_23_scaled) #2133
sample11 <- rXY(n11,estimates1[1:3],c(estimates1[4],estimates1[4]),c(estimates1[5],estimates1[5]),estimates1[6],PLOT=F,onlyXY=FALSE)
sample22 <- rXY(n22,estimates2[1:3],c(estimates2[4],estimates2[4]),c(estimates2[5],estimates2[5]),estimates2[6],PLOT=F,onlyXY=FALSE)
sample33 <- rXY(n33,estimates3[1:3],c(estimates3[4],estimates3[4]),c(estimates3[5],estimates3[5]),estimates3[6],PLOT=F,onlyXY=FALSE)
ssX1 <- sample11$XY[,1]
ssY1 <- sample11$XY[,2]
ssX2 <- sample22$XY[,1]
ssY2 <- sample22$XY[,2]
ssX3 <- sample33$XY[,1]
ssY3 <- sample33$XY[,2]

### pair-plots: real data vs simulated
par(mfrow=c(1,3), mar = c(4,4,4,4))
# First row  
plot(ssX1,ssY1, xlab="AMMERZODEN", ylab="GIERSBERGEN", pch=20, cex=.9)
points(mat_11_15_scaled[1,],mat_11_15_scaled[2,] ,pch=20, cex=.9, col="red")
plot(ssX2,ssY2, xlab="AMMERZODEN", ylab="ZALTBOMMEL", pch=20, cex=.9)
points(mat_11_23_scaled[1,],mat_11_23_scaled[2,] ,pch=20, cex=.9, col="red")
plot(ssX3,ssY3, xlab="GIERSBERGEN", ylab="ZALTBOMMEL", pch=20, cex=.9)
points(mat_15_23_scaled[1,],mat_15_23_scaled[2,] ,pch=20, cex=.9, col="red")

##################################################
##  Chi-measure btw the real and simulated data ##
##################################################


####  chi(u) and chi.bar(u) measures ####
#########################################

chi.u = function(data,U){
  sum(data[,1]>U & data[,2]>U)/sum(data[,2]>U)
}

chi.l = function(data,L){
  sum(data[,1]<L & data[,2]<L)/sum(data[,2]<L)
}


#Empirical Transform to uniform(0,1):
n <- 10^7
sim.XY1 <- rXY(n,estimates1[1:3],c(estimates1[4],estimates1[4]),c(estimates1[5],estimates1[5]),estimates1[6],PLOT=F,onlyXY=FALSE)
sim.XY2 <- rXY(n,estimates2[1:3],c(estimates2[4],estimates2[4]),c(estimates2[5],estimates2[5]),estimates2[6],PLOT=F,onlyXY=FALSE)
sim.XY3 <- rXY(n,estimates3[1:3],c(estimates3[4],estimates3[4]),c(estimates3[5],estimates3[5]),estimates3[6],PLOT=F,onlyXY=FALSE)
X1 <- sim.XY1$XY[,1]
Y1 <- sim.XY1$XY[,2]
X2 <- sim.XY2$XY[,1]
Y2 <- sim.XY2$XY[,2]
X3 <- sim.XY3$XY[,1]
Y3 <- sim.XY3$XY[,2]
sim.11_15 <- cbind(X1,Y1)
sim.11_23 <- cbind(X2,Y2)
sim.15_23 <- cbind(X3,Y3)

X11 <- ecdf(X1)
Y11 <- ecdf(Y1)
X22 <- ecdf(X2)
Y22 <- ecdf(Y2)
X33 <- ecdf(X3)
Y33 <- ecdf(Y3)
sim1.11_15 <- cbind(X11(X1),Y11(Y1))
sim1.11_23 <- cbind(X22(X2),Y22(Y2))
sim1.15_23 <- cbind(X33(X3),Y33(Y3))
data_pair1_11 <- ecdf(mat_11_15_scaled[1,])
data_pair1_15 <- ecdf(mat_11_15_scaled[2,])
data_pair2_11 <- ecdf(mat_11_23_scaled[1,])
data_pair2_23 <- ecdf(mat_11_23_scaled[2,])
data_pair3_15 <- ecdf(mat_15_23_scaled[1,])
data_pair3_23 <- ecdf(mat_15_23_scaled[2,])
ECA_data_11_15 <- cbind(data_pair1_11(mat_11_15_scaled[1,]), data_pair1_15(mat_11_15_scaled[2,]))
ECA_data_11_23 <- cbind(data_pair2_11(mat_11_23_scaled[1,]), data_pair2_23(mat_11_23_scaled[2,]))
ECA_data_15_23 <- cbind(data_pair3_15(mat_15_23_scaled[1,]), data_pair3_23(mat_15_23_scaled[2,]))


U=seq(0,1, by=0.01)
L=seq(0,1, by=0.01)
chi.U1=c()
chi.L1=c()
chi.Udata1=c()
chi.Ldata1=c()
chi.U2=c()
chi.L2=c()
chi.Udata2=c()
chi.Ldata2=c()
chi.U3=c()
chi.L3=c()
chi.Udata3=c()
chi.Ldata3=c()

for (i in 1:length(U)){
  chi.U1[i] <- chi.u(sim1.11_15,U[i])
  chi.Udata1[i] <- chi.u(ECA_data_11_15,U[i])
}
for (i in 1:length(U)){
  chi.U2[i] <- chi.u(sim1.11_23,U[i])
  chi.Udata2[i] <- chi.u(ECA_data_11_23,U[i])
}
for (i in 1:length(U)){
  chi.U3[i] <- chi.u(sim1.15_23,U[i])
  chi.Udata3[i] <- chi.u(ECA_data_15_23,U[i])
}

for (i in 1:length(L)){
  chi.L1[i] <- chi.l(sim1.11_15,L[i])
  chi.Ldata1[i] <- chi.l(ECA_data_11_15,L[i])
}
for (i in 1:length(L)){
  chi.L2[i] <- chi.l(sim1.11_23,L[i])
  chi.Ldata2[i] <- chi.l(ECA_data_11_23,L[i])
}
for (i in 1:length(L)){
  chi.L3[i] <- chi.l(sim1.15_23,L[i])
  chi.Ldata3[i] <- chi.l(ECA_data_15_23,L[i])
}


### plots:
par(mfrow=c(3,2), mar = c(4,4,4,4))
# First row
plot(U,chi.U1, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")))
lines(U,chi.Udata1, col = "red")
plot(L,chi.L1, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")))
lines(L,chi.Ldata1, col = "red")
# Second row
plot(U,chi.U2, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")))
lines(U,chi.Udata2, col = "red")
plot(L,chi.L2, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")))
lines(L,chi.Ldata2, col = "red")
# Third row
plot(U,chi.U3, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")))
lines(U,chi.Udata3, col = "red")
plot(L,chi.L3, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")))
lines(L,chi.Ldata3, col = "red")


####################################################
####   Obtaining a confidence interval of 95%  #####
####################################################

#################################################################
## Bootstrap function: (with pointwise confidence bands)
#################################################################
CBbootstrapping <- function(iter,u,data){
  res <- list()
  boot_chi_u <- c()
  boot_chi_l <- c()
  p <- length(data[,1])
  
  for (i in 1:iter){
    ind <- sample(1:p, p, replace = T)
    data_boot <- data[ind,]
    
    for(j in 1:length(u)){
      boot_chi_u[j] <- chi.u(data_boot, u[j])
      boot_chi_l[j] <- chi.l(data_boot, u[j])
    }
    res[[i]] <- cbind(boot_chi_u, boot_chi_l)
  }
  return(res)
  
}


#########
RNGversion("3.6.0")
set.seed(1001)


U <- seq(0,1, by = 0.01)
L <- seq(0,1, by=0.01)
iter <- 10000
Timestart <- Sys.time()
CBbootstrapping_11_15 <- CBbootstrapping(iter, U, data = ECA_data_11_15)
Timeend2 <- Sys.time() 
TotTime2 <- Timeend2 - Timestart # 3.056346 mins

CBbootstrapping_11_23 <- CBbootstrapping(iter, U, data = ECA_data_11_23)
CBbootstrapping_15_23 <- CBbootstrapping(iter, U, data = ECA_data_15_23)

Boot.CI.u1 <- matrix(NA,ncol=iter,nrow=length(U))
Boot.CI.l1 <- matrix(NA,ncol=iter,nrow=length(U))
Boot.CI.u2 <- matrix(NA,ncol=iter,nrow=length(U))
Boot.CI.l2 <- matrix(NA,ncol=iter,nrow=length(U))
Boot.CI.u3 <- matrix(NA,ncol=iter,nrow=length(U))
Boot.CI.l3 <- matrix(NA,ncol=iter,nrow=length(U))

for (i in 1:iter){
  Boot.CI.u1[,i] <- CBbootstrapping_11_15[[i]][,1]
  Boot.CI.l1[,i] <- CBbootstrapping_11_15[[i]][,2]
}
for (i in 1:iter){
  Boot.CI.u2[,i] <- CBbootstrapping_11_23[[i]][,1]
  Boot.CI.l2[,i] <- CBbootstrapping_11_23[[i]][,2]
}
for (i in 1:iter){
  Boot.CI.u3[,i] <- CBbootstrapping_15_23[[i]][,1]
  Boot.CI.l3[,i] <- CBbootstrapping_15_23[[i]][,2]
}

#Boot.CB.u2 <- data.frame(t(apply(Boot.CI.u2,1,quantile,c(0.025,.975),na.rm=T)))
#Boot.CB.l2 <- data.frame(t(apply(Boot.CI.l2,1,quantile,c(0.025,.975),na.rm=T)))

#################################################################
##### Using envelope to compute point-wise confidence bands: ####
#################################################################
library(boot)
Timestart <- Sys.time()
Boot.CI.u1 <- Boot.CI.u1[-101,]
Boot.CI.l1 <- Boot.CI.l1[-1,]
PW.CB.u1 <- envelope(mat = t(Boot.CI.u1), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.u1)))
PW.CB.l1 <- envelope(mat = t(Boot.CI.l1), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.l1)))
Timeend2 <- Sys.time() 
TotTime2 <- Timeend2 - Timestart

Timestart <- Sys.time()
Boot.CI.u2 <- Boot.CI.u2[-101,]
Boot.CI.l2 <- Boot.CI.l2[-1,]
PW.CB.u2 <- envelope(mat = t(Boot.CI.u2), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.u2)))
PW.CB.l2 <- envelope(mat = t(Boot.CI.l2), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.l2)))
Timeend2 <- Sys.time() 
TotTime2 <- Timeend2 - Timestart #1.527398 secs

Boot.CI.u3 <- Boot.CI.u3[-101,]
Boot.CI.l3 <- Boot.CI.l3[-1,]
PW.CB.u3 <- envelope(mat = t(Boot.CI.u3), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.u3)))
PW.CB.l3 <- envelope(mat = t(Boot.CI.l3), level = c(0.95,0.95), index = 1:ncol(t(Boot.CI.l3)))

# Xi-measures plots using envelope: 
# first row: Data 11_15
par(mfrow=c(3,2), mar = c(4,4,4,4))
plot(U,chi.U1, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")), col = "red", asp=1)
lines(U,chi.Udata1)
lines(U[-101], PW.CB.u1$point[1, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u1$point[2, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u1$overall[1, ], lty = 1, col = "gray")
lines(U[-101], PW.CB.u1$overall[2, ], lty = 1, col = "gray")


plot(L,chi.L1, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")), col = "red", asp=1)
lines(L,chi.Ldata1)
lines(L[-1], PW.CB.l1$point[1, ], lty = 4, col = "gray")
lines(L[-1], PW.CB.l1$point[2, ], lty = 4, col = "gray")
lines(L[-1], PW.CB.l1$overall[1, ], lty = 1, col = "gray")
lines(L[-1], PW.CB.l1$overall[2, ], lty = 1, col = "gray")

# second row: Data 11_23
plot(U,chi.U2, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")), col = "red", asp=1)
lines(U,chi.Udata2)
lines(U[-101], PW.CB.u2$point[1, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u2$point[2, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u2$overall[1, ], lty = 1, col = "gray")
lines(U[-101], PW.CB.u2$overall[2, ], lty = 1, col = "gray")


plot(L,chi.L2, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")), col = "red", asp=1)
lines(L,chi.Ldata2)
lines(L[c(-1,-2)], PW.CB.l2$point[1, ], lty = 4, col = "gray")
lines(L[c(-1,-2)], PW.CB.l2$point[2, ], lty = 4, col = "gray")
lines(L[c(-1,-2)], PW.CB.l2$overall[1, ], lty = 1, col = "gray")
lines(L[c(-1,-2)], PW.CB.l2$overall[2, ], lty = 1, col = "gray")

# third row: Data 15_23 
plot(U,chi.U3, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold u", ylab=expression(paste(chi,"(u)")), col = "red", asp=1)
lines(U,chi.Udata3)
lines(U[-101], PW.CB.u3$point[1, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u3$point[2, ], lty = 4, col = "gray")
lines(U[-101], PW.CB.u3$overall[1, ], lty = 1, col = "gray")
lines(U[-101], PW.CB.u3$overall[2, ], lty = 1, col = "gray")


plot(L,chi.L3, type="l",xlim=c(0,1),ylim=c(0,1), xlab="Threshold l", ylab=expression(paste(chi,"(l)")), col = "red", asp=1)
lines(L,chi.Ldata3)
lines(L[c(-1,-2)], PW.CB.l3$point[1, ], lty = 4, col = "gray")
lines(L[c(-1,-2)], PW.CB.l3$point[2, ], lty = 4, col = "gray")
lines(L[c(-1,-2)], PW.CB.l3$overall[1, ], lty = 1, col = "gray")
lines(L[c(-1,-2)], PW.CB.l3$overall[2, ], lty = 1, col = "gray")


#################################################################
##### Compute point-wise confidence bands for the QQ plots:  ####
#################################################################

QQbootstrapping <- function(iter, m, data){
  res <- list()
  boot1 <- c()
  boot2 <- c()
  boot3 <- c()
  p <- length(data[,1])
  
  for (i in 1:iter){
    ind <- sample(1:p, p, replace = T)
    data_boot <- data[ind,]
    
    for(j in 1:m){
      boot1[j] <- as.numeric(quantile(data_boot[,1], j/(m+1), na.rm = T))
      boot2[j] <- as.numeric(quantile(data_boot[,2], j/(m+1), na.rm = T))
      boot3[j] <- as.numeric(quantile(data_boot[,3], j/(m+1), na.rm = T))
    }
    res[[i]] <- cbind(boot1, boot2, boot3)
  }
  return(res)
  
}

RNGversion("3.6.0")
set.seed(1001)


ssXY1 <- ssX1 + ssY1
ssXY2 <- ssX2 + ssY2
ssXY3 <- ssX3 + ssY3

iter <- 10000
Timestart <- Sys.time()
QQbootstrapping1 <- QQbootstrapping(iter, n11, data = cbind(ssX1,ssY1,ssXY1))
Timeend2 <- Sys.time()
TotTime2 <- Timeend2 - Timestart #4.5 hours

Timestart <- Sys.time()
QQbootstrapping2 <- QQbootstrapping(iter, n22, data = cbind(ssX2,ssY2,ssXY2))
Timeend2 <- Sys.time()
TotTime2 <- Timeend2 - Timestart  #4.2 hours

Timestart2 <- Sys.time()
QQbootstrapping3 <- QQbootstrapping(iter, n33, data = cbind(ssX3,ssY3,ssXY3))
Timeend22 <- Sys.time()
TotTime22 <- Timeend22 - Timestart2 #3.5 hours

QQBoot.CI.X1 <- matrix(NA,ncol=iter, nrow=n11)
QQBoot.CI.Y1 <- matrix(NA, ncol=iter, nrow=n11)
QQBoot.CI.XY1 <- matrix(NA, ncol=iter, nrow=n11)
QQBoot.CI.X2 <- matrix(NA,ncol=iter, nrow=n22)
QQBoot.CI.Y2 <- matrix(NA, ncol=iter, nrow=n22)
QQBoot.CI.XY2 <- matrix(NA, ncol=iter, nrow=n22)
QQBoot.CI.X3 <- matrix(NA,ncol=iter, nrow=n33)
QQBoot.CI.Y3 <- matrix(NA, ncol=iter, nrow=n33)
QQBoot.CI.XY3 <- matrix(NA, ncol=iter, nrow=n33)

for (i in 1:iter){
  QQBoot.CI.X1[,i] <- QQbootstrapping1[[i]][,1]
  QQBoot.CI.Y1[,i] <- QQbootstrapping1[[i]][,2]
  QQBoot.CI.XY1[,i] <- QQbootstrapping1[[i]][,3]
  QQBoot.CI.X2[,i] <- QQbootstrapping2[[i]][,1]
  QQBoot.CI.Y2[,i] <- QQbootstrapping2[[i]][,2]
  QQBoot.CI.XY2[,i] <- QQbootstrapping2[[i]][,3]
  QQBoot.CI.X3[,i] <- QQbootstrapping3[[i]][,1]
  QQBoot.CI.Y3[,i] <- QQbootstrapping3[[i]][,2]
  QQBoot.CI.XY3[,i] <- QQbootstrapping3[[i]][,3]
}

save(
  QQBoot.CI.X2,
  QQBoot.CI.Y2,
  QQBoot.CI.XY2,
  file = "QQ_bootstrap_pair11_23.RData"
)

save(
  QQBoot.CI.X1,
  QQBoot.CI.Y1,
  QQBoot.CI.XY1,
  file = "QQ_bootstrap_pair11_15.RData"
)

save(
  QQBoot.CI.X3,
  QQBoot.CI.Y3,
  QQBoot.CI.XY3,
  file = "QQ_bootstrap_pair15_23.RData"
)

# for pair data 2
Boot.CB.x2 <- data.frame(t(apply(QQBoot.CI.X2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.x2.mean <- data.frame(apply(QQBoot.CI.X2,1,mean,na.rm=T))
Boot.CB.x2.med <- data.frame(apply(QQBoot.CI.X2,1,quantile,0.5,na.rm=T))
Boot.CB.y2 <- data.frame(t(apply(QQBoot.CI.Y2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.y2.mean <- data.frame(apply(QQBoot.CI.Y2,1,mean,na.rm=T))
Boot.CB.y2.med <- data.frame(apply(QQBoot.CI.Y2,1,quantile,0.5,na.rm=T))
Boot.CB.xy2 <- data.frame(t(apply(QQBoot.CI.XY2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.xy2.mean <- data.frame(apply(QQBoot.CI.XY2,1,mean,na.rm=T))
Boot.CB.xy2.med <- data.frame(apply(QQBoot.CI.XY2,1,quantile,0.5,na.rm=T))

# then in any future session:
#load("QQ_bootstrap_pair11_23.RData")
#load("QQ_bootstrap_pair11_15.RData")
#load("QQ_bootstrap_pair15_23.RData")


# for pair data 1
Boot.CB.x1.med <- data.frame(apply(QQBoot.CI.X1,1,quantile,0.5,na.rm=T))
Boot.CB.y1.med <- data.frame(apply(QQBoot.CI.Y1,1,quantile,0.5,na.rm=T))
Boot.CB.xy1.med <- data.frame(apply(QQBoot.CI.XY1,1,quantile,0.5,na.rm=T))
Boot.CB.x1.mean <- data.frame(apply(QQBoot.CI.X1,1,mean,na.rm=T))
Boot.CB.y1.mean <- data.frame(apply(QQBoot.CI.Y1,1,mean,na.rm=T))
Boot.CB.xy1.mean <- data.frame(apply(QQBoot.CI.XY1,1,mean,na.rm=T))

# for pair data 2
Boot.CB.x2 <- data.frame(t(apply(QQBoot.CI.X2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.x2.mean <- data.frame(apply(QQBoot.CI.X2,1,mean,na.rm=T))
Boot.CB.x2.med <- data.frame(apply(QQBoot.CI.X2,1,quantile,0.5,na.rm=T))
Boot.CB.y2 <- data.frame(t(apply(QQBoot.CI.Y2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.y2.mean <- data.frame(apply(QQBoot.CI.Y2,1,mean,na.rm=T))
Boot.CB.y2.med <- data.frame(apply(QQBoot.CI.Y2,1,quantile,0.5,na.rm=T))
Boot.CB.xy2 <- data.frame(t(apply(QQBoot.CI.XY2,1,quantile,c(0.025,.975),na.rm=T)))
Boot.CB.xy2.mean <- data.frame(apply(QQBoot.CI.XY2,1,mean,na.rm=T))
Boot.CB.xy2.med <- data.frame(apply(QQBoot.CI.XY2,1,quantile,0.5,na.rm=T))

# for pair data 3
Boot.CB.x3.med <- data.frame(apply(QQBoot.CI.X3,1,quantile,0.5,na.rm=T))
Boot.CB.y3.med <- data.frame(apply(QQBoot.CI.Y3,1,quantile,0.5,na.rm=T))
Boot.CB.xy3.med <- data.frame(apply(QQBoot.CI.XY3,1,quantile,0.5,na.rm=T))
Boot.CB.x3.mean <- data.frame(apply(QQBoot.CI.X3,1,mean,na.rm=T))
Boot.CB.y3.mean <- data.frame(apply(QQBoot.CI.Y3,1,mean,na.rm=T))
Boot.CB.xy3.mean <- data.frame(apply(QQBoot.CI.XY3,1,mean,na.rm=T))

### Using envolpe package to compute overall error bound + pairwise confidance bounds
PW.CB.x1 <- envelope(mat = t(QQBoot.CI.X1), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.X1)))
PW.CB.y1 <- envelope(mat = t(QQBoot.CI.Y1), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.Y1)))
PW.CB.xy1 <- envelope(mat = t(QQBoot.CI.XY1), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.XY1)))

PW.CB.x2 <- envelope(mat = t(QQBoot.CI.X2), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.X2)))
PW.CB.y2 <- envelope(mat = t(QQBoot.CI.Y2), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.Y2)))
PW.CB.xy2 <- envelope(mat = t(QQBoot.CI.XY2), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.XY2)))

PW.CB.x3 <- envelope(mat = t(QQBoot.CI.X3), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.X3)))
PW.CB.y3 <- envelope(mat = t(QQBoot.CI.Y3), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.Y3)))
PW.CB.xy3 <- envelope(mat = t(QQBoot.CI.XY3), level = c(0.95,0.95), index = 1:ncol(t(QQBoot.CI.XY3)))

###################### for data pair 1 ######################   
plot(sort(mat_11_15_scaled[1,]), Boot.CB.x1.mean[,1], xlab = "simulated X", ylab = "AMMERZODEN",  main = "")
abline(0, 1, col = "lightgrey")
lines(sort(mat_11_15_scaled[1,]), PW.CB.x1$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[1,]), PW.CB.x1$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[1,]), PW.CB.x1$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_15_scaled[1,]), PW.CB.x1$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_11_15_scaled[2,]), Boot.CB.y1.mean[,1], xlab = "simulated Y", ylab = "GIERSBERGEN",  main = "")
abline(0, 1, col = "lightgrey")
lines(sort(mat_11_15_scaled[2,]), PW.CB.y1$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[2,]), PW.CB.y1$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[2,]), PW.CB.y1$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_15_scaled[2,]), PW.CB.y1$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_11_15_scaled[1,]+mat_11_15_scaled[2,]), Boot.CB.xy1.mean[,1], xlab = "simulated X+Y", ylab = "AMMERZODEN + GIERSBERGEN",  main = "")
abline(0, 1, col = "lightgrey") 
lines(sort(mat_11_15_scaled[1,]+mat_11_15_scaled[2,]), PW.CB.xy1$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[1,]+mat_11_15_scaled[2,]), PW.CB.xy1$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_15_scaled[1,]+mat_11_15_scaled[2,]), PW.CB.xy1$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_15_scaled[1,]+mat_11_15_scaled[2,]), PW.CB.xy1$overall[2, ], lty = 1, col = "gray")


###################### for data pair 2 ######################
par(mfrow=c(1,3), xpd   = FALSE )

plot(sort(mat_11_23_scaled[1,]), Boot.CB.x2.mean[,1], xlab = "simulated X", ylab = "AMMERZODEN",  main = "")
abline(0, 1, col = "lightgrey")
lines(sort(mat_11_23_scaled[1,]), PW.CB.x2$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[1,]), PW.CB.x2$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[1,]), PW.CB.x2$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_23_scaled[1,]), PW.CB.x2$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_11_23_scaled[2,]), Boot.CB.y2.mean[,1], xlab = "simulated Y", ylab = "ZALTBOMMEL",  main = "")
abline(0, 1, col = "lightgrey") 
lines(sort(mat_11_23_scaled[2,]), PW.CB.y2$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[2,]), PW.CB.y2$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[2,]), PW.CB.y2$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_23_scaled[2,]), PW.CB.y2$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), Boot.CB.xy2.mean[,1], xlab = "simulated X+Y", ylab = "AMMERZODEN + ZALTBOMMEL",  main = "")
abline(0, 1, col = "lightgrey") 
lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy2$point[1, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy2$point[2, ], lty = 4, col = "gray")
lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy2$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy2$overall[2, ], lty = 1, col = "gray")

###################### for data pair 3 ######################
plot(sort(mat_15_23_scaled[1,]), Boot.CB.x3.mean[,1], xlab = "simulated X", ylab = "GIERSBERGEN",  main = "")
abline(0, 1, col = "lightgrey")
lines(sort(mat_15_23_scaled[1,]), PW.CB.x3$point[1, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[1,]), PW.CB.x3$point[2, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[1,]), PW.CB.x3$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_15_23_scaled[1,]), PW.CB.x3$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_15_23_scaled[2,]), Boot.CB.y3.mean[,1], xlab = "simulated Y", ylab = "ZALTBOMMEL",  main = "")
abline(0, 1, col = "lightgrey")
lines(sort(mat_15_23_scaled[2,]), PW.CB.y3$point[1, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[2,]), PW.CB.y3$point[2, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[2,]), PW.CB.y3$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_15_23_scaled[2,]), PW.CB.y3$overall[2, ], lty = 1, col = "gray")

plot(sort(mat_15_23_scaled[1,]+mat_15_23_scaled[2,]), Boot.CB.xy3.mean[,1], xlab = "simulated X+Y", ylab = "GIERSBERGEN + ZALTBOMMEL",  main = "")
abline(0, 1, col = "lightgrey") 
lines(sort(mat_15_23_scaled[1,]+mat_15_23_scaled[2,]), PW.CB.xy3$point[1, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[1,]+mat_15_23_scaled[2,]), PW.CB.xy3$point[2, ], lty = 4, col = "gray")
lines(sort(mat_15_23_scaled[1,]+mat_15_23_scaled[2,]), PW.CB.xy3$overall[1, ], lty = 1, col = "gray")
lines(sort(mat_15_23_scaled[1,]+mat_15_23_scaled[2,]), PW.CB.xy3$overall[2, ], lty = 1, col = "gray")


