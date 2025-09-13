library("MASS")
library("dplyr")
library("parallel")
library("NeuralEstimators")
library("ggplot2")
source("R/simulations/ModelSimulator.R")

NaiveEstimator <- function(Y,thresh,start,param.bounds){
  threshL <- thresh[1]
  threshU <- thresh[2]
  
  start.kappa <- start[1]
  start.sigma <- start[2]
  start.xi <- start[3]
  start.thL <- start[4]
  start.thU <- start[5]
  start.thw <- start[6]
  
  ## ML Estimation of kappa,sigma,xi
  Y.norm <- apply(Y,1,sum)
  
  neg.log.likelihood.eGPD <- function(param,R){
    kappa <- param[1]
    sigma <- param[2]
    xi <- param[3]
    if(kappa>0 & sigma>0 & xi>0 & !any(R<=0)){
      vals <- -log(kappa) - (kappa-1)*log(1-(1+xi*R/sigma)^(-1/xi)) + (1/xi+1)*log(1+xi*R/sigma) + log(sigma)
      res <- sum(vals[!is.na(vals) & vals!=Inf])
      return(res)
    } else{
      return( 10^9 )
    }
  }
  eGPD.fit <- tryCatch({optim(
    par=c(start.kappa,start.sigma,start.xi),
    fn=neg.log.likelihood.eGPD,
    R=Y.norm,
    method="L-BFGS-B",
    lower=param.bounds[1:3,1],
    upper=param.bounds[1:3,2],
    hessian=FALSE,
    control=list(maxit=1000)
    )$par
  }, error = function(e) {
    warning(paste("Optimization failed. Error:", e$message))
    return(start)
  })
  
  kappa.est <- eGPD.fit[1]
  sigma.est <- eGPD.fit[2]
  xi.est <- eGPD.fit[3]  
  
  neg.log.likelihood.beta <- function(param,data){
    th <- param[1]
    if(th>0){
      res <- -sum(dbeta(data,th,th,log=TRUE)) 
      if (is.infinite(res)) {
        return(10^9)
      } else {
        return(res) 
      }
    } else{
      return( 10^9 )
    }
  }
  
  angleL <- Y[Y.norm<quantile(Y.norm,threshL) & Y.norm>0,]/Y.norm[Y.norm<quantile(Y.norm,threshL) & Y.norm>0]
  betaL.fit <- optim(par=start.thL,fn=neg.log.likelihood.beta,data=angleL[,1],method="L-BFGS-B",lower=param.bounds[4,1],upper=param.bounds[4,2],hessian=FALSE,control=list(maxit=1000))
  thL.est <- betaL.fit$par
  
  angleU <- Y[Y.norm>quantile(Y.norm,threshU),]/Y.norm[Y.norm>quantile(Y.norm,threshU)]
  betaU.fit <- optim(par=start.thU,fn=neg.log.likelihood.beta,data=angleU[,1],method="L-BFGS-B",lower=param.bounds[5,1],upper=param.bounds[5,2],hessian=FALSE,control=list(maxit=1000))
  thU.est <- betaU.fit$par
  
  LS.omega <- function(param,data,param.fixed,Nsamples=10^6){
    if(param>0 & param<0.5){
      set.seed(8349862)
      samples <- rsim(Nsamples,c(param.fixed,param))
      res <- sum((cov(samples$Y)-cov(data))^2)
      return(res)
    } else{
      return( 10^9 )
    }
  }
  omega.fit <- optim(par=start.thw,fn=LS.omega,data=Y,param.fixed=c(kappa.est,sigma.est,xi.est,thL.est,thU.est),Nsamples=10^6,method="L-BFGS-B",lower=param.bounds[6,1],upper=param.bounds[6,2],hessian=FALSE,control=list(maxit=1000))
  thw.est <- omega.fit$par
  
  return(c(kappa.est,sigma.est,xi.est,thL.est,thU.est,thw.est))
}


int_path <- "intermediates"
dir.create(int_path, recursive = TRUE, showWarnings = FALSE)
theta_test <- readRDS("intermediates/theta_test.rds")
Z_test <- readRDS("intermediates/Z_test.rds")

thresh <- c(0.05,0.95)
param.bounds <- matrix(c(0.1,0.1,0,0.1,0.1,0.01,
                         20,3,0.5,20,20,0.49),
                       nrow=6)

n.theta_test <- ncol(theta_test)

theta.hat.Naive <- matrix(nrow=6,ncol=n.theta_test)
rownames(theta.hat.Naive) <- c("kappa", "sigma", "xi", "thL", "thU", "thw")
for(i in 1:n.theta_test){ 
  print(i)
  if(i==124 | i==796 | i==817){ ## for these cases, the optimization fails if we use the "prior mean" as initial values...
    start <- theta_test[,i] ## Initial values are true values
  } else{
    start <- apply(param.bounds,1,diff)/2 ## Initial values are the "prior mean"
  }
  theta.hat.Naive[,i] <- NaiveEstimator(Y=t(Z_test[[i]]),thresh=thresh,start=start,param.bounds=param.bounds)
  cat("truth:   ", round(theta_test[,i], 3), "\n")
  cat("estimate:", round(theta.hat.Naive[,i], 3), "\n")
}
saveRDS(theta.hat.Naive, file = file.path(int_path, "theta_test_naiveestimates.rds"))

df.naive <- data.frame(
  estimator = rep("Naive",length(theta_test)), 
  parameter = rep(rownames(theta_test),ncol(theta_test)), 
  truth = as.vector(theta_test), 
  estimate = as.vector(theta.hat.Naive)
)

## Save figure into intermediates/Figures/Simulation
fig_path <- file.path("intermediates", "Figures", "Simulation")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

ggsave(
  plot = plotestimates(df.naive), 
  file = "Naive_assessment.pdf", 
  path = fig_path, device = "pdf", width = 18, height = 3
)


## Inference time:
inference_time <- system.time({
  i = 1
  start <- theta_test[,i] ## Initial values are true values
  NaiveEstimator(Y=t(Z_test[[1]]),thresh=thresh,start=start,param.bounds=param.bounds)
})
writeLines(as.character(inference_time["elapsed"]), file.path(int_path, "inference_time_Naive.txt"))
