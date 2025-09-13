###############################################################
## ECA_data_diagnostic_plots.R – Diagnostics for Pair 11–23 ##
###############################################################

library("NeuralEstimators")
library("JuliaConnectoR")
library("boot")        # for envelope()
source("R/application/MEGPD_Model.R")
#Sys.setenv(JULIA_BINDIR = "/Users/alotainm/.julia/juliaup/julia-1.11.5+0.x64.apple.darwin14/bin")
#Sys.setenv("JULIACONNECTOR_JULIAOPTS" = "--project=.")
#juliaEval('using NeuralEstimators, Flux')
#source("R/simulations/Architecture.R")
# Functions from Julia
#sampleposterior <- juliaFun("sampleposterior")
#logdensity <- juliaFun("logdensity")


###############################################################
## 1. Load preprocessed rainfall data (scaled pairs)
###############################################################
df_11_23_scaled <- read.csv("R/data/pair_11_23_matrix_scaled.csv",
                            stringsAsFactors = FALSE)
mat_11_23_scaled <- as.matrix(df_11_23_scaled[,-1])  # drop date column if present
dim(mat_11_23_scaled)  # (2 × 2166)

###############################################################
## 2. Load trained estimator & sample posterior
###############################################################
set.seed(123)
#loadstate(NPE, file.path("intermediates", "NPE.bson"))

# Variance stabilizing transformation
#signed_log <- function(x) (sign(x) * log1p(abs(x))) - 1
#mat_11_23 <- t(apply(mat_11_23_scaled, 1, signed_log))

# Posterior samples
#samples_11_23 <- sampleposterior(NPE, mat_11_23)[[1]]
#samples_11_23 <- exp(samples_11_23)  # back-transform
#estimates <- apply(samples_11_23, 1, median)
#estimates_CI <- apply(samples_11_23, 1, quantile, c(0.025, 0.975))
estimates <- c(1.1188870,  1.3769174,  0.1953660,  4.0885365, 20.9920100,  0.2361101) 

###############################################################
## 3. Simulate data under estimated parameters
###############################################################
n2 <- ncol(mat_11_23_scaled) # 2166 obs
sample22 <- rXY(n2,
                estimates[1:3],
                c(estimates[4],estimates[4]),
                c(estimates[5],estimates[5]),
                estimates[6],
                PLOT=FALSE, onlyXY=FALSE)

ssX2 <- sample22$XY[,1]
ssY2 <- sample22$XY[,2]

### pair-plots: real data vs simulated
png("intermediates/Figures/Application/scatterplot_pair11_23.png", width=7, height=5, units="in", res=300)
par(mfrow=c(1,1))
plot(ssX2,ssY2, xlab="AMMERZODEN", ylab="ZALTBOMMEL", pch=20, cex=.9)
points(mat_11_23_scaled[1,],mat_11_23_scaled[2,] ,pch=20, cex=.9, col="red")
dev.off()
print("scatter plot done")
###############################################################
## 4. χ(u) and χ(l) measures with bootstrap CIs
###############################################################
chi.u <- function(data,U) sum(data[,1]>U & data[,2]>U)/sum(data[,2]>U)
chi.l <- function(data,L) sum(data[,1]<L & data[,2]<L)/sum(data[,2]<L)

# Transform sim & data to uniforms
sim_big <- rXY(1e6, estimates[1:3],
               c(estimates[4],estimates[4]),
               c(estimates[5],estimates[5]),
               estimates[6], PLOT=FALSE, onlyXY=FALSE)

X2 <- sim_big$XY[,1]; Y2 <- sim_big$XY[,2]
sim1.11_23 <- cbind(ecdf(X2)(X2), ecdf(Y2)(Y2))

data_pair2_11 <- ecdf(mat_11_23_scaled[1,])
data_pair2_23 <- ecdf(mat_11_23_scaled[2,])
ECA_data_11_23 <- cbind(data_pair2_11(mat_11_23_scaled[1,]),
                        data_pair2_23(mat_11_23_scaled[2,]))

# Compute chi(u), chi(l)
U <- seq(0,1, by=0.01)
L <- seq(0,1, by=0.01)
iter <- 100   # reduce if too slow
chi.U <- sapply(U, function(u) chi.u(sim1.11_23,u))
chi.Udata <- sapply(U, function(u) chi.u(ECA_data_11_23,u))
chi.L <- sapply(L, function(l) chi.l(sim1.11_23,l))
chi.Ldata <- sapply(L, function(l) chi.l(ECA_data_11_23,l))

# Bootstrap envelopes
CBbootstrapping <- function(iter,u,data){
  res <- vector("list", iter)
  for (i in 1:iter){
    ind <- sample(1:nrow(data), nrow(data), replace=TRUE)
    data_boot <- data[ind,]
    boot_chi_u <- sapply(u, function(val) chi.u(data_boot,val))
    boot_chi_l <- sapply(u, function(val) chi.l(data_boot,val))
    res[[i]] <- cbind(boot_chi_u, boot_chi_l)
  }
  res
}

RNGversion("3.6.0") 
set.seed(1001)
iter <- 100   # reduce if too slow
boot_res <- CBbootstrapping(iter, U, data=ECA_data_11_23)

Boot.CI.u <- sapply(boot_res, function(x) x[,1])
Boot.CI.l <- sapply(boot_res, function(x) x[,2])

Boot.CI.u <- Boot.CI.u[-101, , drop=FALSE] 
Boot.CI.l <- Boot.CI.l[-c(1,2), , drop = FALSE]
                    
PW.CB.u <- envelope(mat=t(Boot.CI.u), level=c(0.95,0.95), index=1:length(U))
PW.CB.l <- envelope(mat=t(Boot.CI.l), level=c(0.95,0.95), index=1:length(L))

# Plot χ(u), χ(l)
png("intermediates/Figures/Application/chi_measures_pair11_23.png", width=7, height=5, units="in", res=300)
par(mfrow=c(1,2))
plot(U, chi.U, type="l", col="red", ylim=c(0,1),
     xlab="Threshold u", ylab=expression(chi(u)))
lines(U, chi.Udata)
lines(U, PW.CB.u$point[1,], lty=2, col="gray")
lines(U, PW.CB.u$point[2,], lty=2, col="gray")

plot(L, chi.L, type="l", col="red", ylim=c(0,1),
     xlab="Threshold l", ylab=expression(chi(l)))
lines(L, chi.Ldata)
lines(L, PW.CB.l$point[1,], lty=2, col="gray")
lines(L, PW.CB.l$point[2,], lty=2, col="gray")
dev.off()

print("chi plot done")

###############################################################
## 5. QQ plot with bootstrap CI
###############################################################
QQbootstrapping <- function(iter, m, data){
  res <- vector("list", iter)
  for (i in 1:iter){
    ind <- sample(1:nrow(data), nrow(data), replace=TRUE)
    data_boot <- data[ind,]
    boot1 <- quantile(data_boot[,1], probs=seq(1,m)/(m+1))
    boot2 <- quantile(data_boot[,2], probs=seq(1,m)/(m+1))
    boot3 <- quantile(rowSums(data_boot), probs=seq(1,m)/(m+1))
    res[[i]] <- cbind(boot1, boot2, boot3)
  }
  res
}

n2 <- ncol(mat_11_23_scaled)
QQboot <- QQbootstrapping(iter, n2, data=cbind(ssX2,ssY2,ssX2+ssY2))

QQBoot.CI.X <- sapply(QQboot, function(x) x[,1])
QQBoot.CI.Y <- sapply(QQboot, function(x) x[,2])
QQBoot.CI.XY <- sapply(QQboot, function(x) x[,3])

PW.CB.x  <- envelope(mat = t(QQBoot.CI.X),  level=c(0.95,0.95), index=1:ncol(t(QQBoot.CI.X)))
PW.CB.y  <- envelope(mat = t(QQBoot.CI.Y),  level=c(0.95,0.95), index=1:ncol(t(QQBoot.CI.Y)))
PW.CB.xy <- envelope(mat = t(QQBoot.CI.XY), level=c(0.95,0.95), index=1:ncol(t(QQBoot.CI.XY)))

# QQ plots
png("intermediates/Figures/Application/qq_plots_pair11_23.png", width=9, height=4, units="in", res=300)
par(mfrow=c(1,3))
plot(sort(mat_11_23_scaled[1,]), rowMeans(QQBoot.CI.X),
     xlab="Simulated X", ylab="Observed (Ammerzoden)")
abline(0,1,col="gray"); lines(sort(mat_11_23_scaled[1,]), PW.CB.x$point[1,], lty=2, col="gray")
lines(sort(mat_11_23_scaled[1,]), PW.CB.x$point[2,], lty=2, col="gray")

plot(sort(mat_11_23_scaled[2,]), rowMeans(QQBoot.CI.Y),
     xlab="Simulated Y", ylab="Observed (Zaltbommel)")
abline(0,1,col="gray"); lines(sort(mat_11_23_scaled[2,]), PW.CB.y$point[1,], lty=2, col="gray")
lines(sort(mat_11_23_scaled[2,]), PW.CB.y$point[2,], lty=2, col="gray")

plot(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), rowMeans(QQBoot.CI.XY),
     xlab="Simulated X+Y", ylab="Observed (Sum)")
abline(0,1,col="gray"); lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy$point[1,], lty=2, col="gray")
lines(sort(mat_11_23_scaled[1,]+mat_11_23_scaled[2,]), PW.CB.xy$point[2,], lty=2, col="gray")
dev.off()

message("Diagnostics for pair 11–23 completed. Plots saved in intermediates/Figures/Application/")

