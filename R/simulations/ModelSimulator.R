rsim <- function(n,param){
  kappa <- param[1]
  sigma <- param[2]
  xi <- param[3]
  thL <- param[4]
  thU <- param[5]
  thw <- param[6]
  
  R.unif <- runif(n)
  R <- sigma*((1-R.unif^(1/kappa))^(-xi)-1)/xi
  V1 <- rbeta(n,thL,thL)
  L1 <- 1/V1
  L2 <- 1/(1-V1)
  L <- cbind(L1,L2)
  U1 <- rbeta(n,thU,thU)
  U2 <- 1-U1
  U <- cbind(U1,U2)
  # wR <- pbeta(R.unif,thw,thw)
  wR <- pbeta((R.unif-thw)/(1-2*thw),3,3)
  
  Y1 <- R*((1-wR)*L1^(-1)+wR*U1)
  Y2 <- R*((1-wR)*L2^(-1)+wR*U2)
  Y <- cbind(Y1,Y2)
  
  return(list(Y=Y,R.unif=R.unif,R=R,L=L,U=U))
}