rm(list=ls())

pEGPD <- function(x,kappa,sigma,xi, LOG=FALSE){
  if(LOG==FALSE){
    res <- (1-(1+xi*x/sigma)^(-1/xi))^kappa
  }else{
    res = kappa*log(1-(1+xi*x/sigma)^(-1/xi))
  }
  return( res )
}

dEGPD <- function(x,kappa,sigma,xi,LOG=FALSE){
  if(LOG==FALSE){
    res <- (kappa/sigma)*(1-(1+xi*x/sigma)^(-1/xi))^(kappa-1)*(1+xi*x/sigma)^(-1/xi-1)
  } else{
    res <- log(kappa)-log(sigma)+(kappa-1)*log(1-(1+xi*x/sigma)^(-1/xi))+(-1/xi-1)*log(1+xi*x/sigma)
  }
  return( res )
}

qEGPD <- function(p,kappa,sigma,xi){
  res <- sigma*((1-p^(1/kappa))^(-xi)-1)/xi
  return( res )
}

rEGPD <- function(n,kappa,sigma,xi){
  res <- qEGPD(runif(n),kappa,sigma,xi)
  return( res )
}

pR <- function(r,parR, LOG=FALSE){
  kappa <- parR[1]
  sigma <- parR[2]
  xi <- parR[3]
  
  res <- pEGPD(r,kappa,sigma,xi,LOG=LOG)
  return( res )
}

dR <- function(r,parR,LOG=FALSE){
  kappa <- parR[1]
  sigma <- parR[2]
  xi <- parR[3]
  
  res <- dEGPD(r,kappa,sigma,xi,LOG=LOG)
  return( res )
}

qR <- function(p,parR){
  kappa <- parR[1]
  sigma <- parR[2]
  xi <- parR[3]
  
  res <- qEGPD(p,kappa,sigma,xi)
  return( res )
}

rR <- function(n,parR){
  kappa <- parR[1]
  sigma <- parR[2]
  xi <- parR[3]
  
  res <- rEGPD(n,kappa,sigma,xi)
  return( res )
}

pWl <- function(w,parWl){
  alpha <- parWl[1]
  beta <- parWl[2]
  
  res <- pbeta(w,shape1=alpha,shape2=beta)
  return( res )
}

dWl <- function(w,parWl,LOG=FALSE){
  alpha <- parWl[1]
  beta <- parWl[2]
  
  if(LOG==FALSE){
    res <- dbeta(w,shape1=alpha,shape2=beta)
  } else{
    res <- dbeta(w,shape1=alpha,shape2=beta,log=TRUE) 
  }
  return( res )
}

qWl <- function(p,parWl){
  alpha <- parWl[1]
  beta <- parWl[2]
  
  res <- qbeta(p,shape1=alpha,shape2=beta)
  return( res )
}

rWl <- function(n,parWl){
  alpha <- parWl[1]
  beta <- parWl[2]
  
  res <- rbeta(n,shape1=alpha,shape2=beta)
  return( res )
}

pWu <- function(w,parWu){
  alpha <- parWu[1]
  beta <- parWu[2]
  
  res <- pbeta(w,shape1=alpha,shape2=beta)
  return( res )
}

dWu <- function(w,parWu,LOG=FALSE){
  alpha <- parWu[1]
  beta <- parWu[2]
  
  res <- dbeta(w,shape1=alpha,shape2=beta,log=LOG)
  return( res )
}

qWu <- function(p,parWu){
  alpha <- parWu[1]
  beta <- parWu[2]
  
  res <- qbeta(p,shape1=alpha,shape2=beta)
  return( res )
}

rWu <- function(n,parWu){
  alpha <- parWu[1]
  beta <- parWu[2]
  
  res <- rbeta(n,shape1=alpha,shape2=beta)
  return( res )
}

weight <- function(cdfR,parWeight){
  alpha <- parWeight[1]
  
  
  res <- pbeta(cdfR,shape1=alpha,shape2=alpha)
  return( res )
}

dL <- function(w,parWu){
  alpha <- parWu[1]
  beta <- parWu[2]
  
  res <- dbeta(1/w,shape1=alpha,shape2=beta)*(1/w^2)
  return( res )
}


rXY <- function(n,parR,parWl,parWu,parWeight,PLOT=TRUE,onlyXY=FALSE){
  R <- rR(n,parR)
  Wl <- rWl(n,parWl)
  Wu <- rWu(n,parWu)
  cdfR <- pR(R,parR) #on uniform scale (defined on [0,1])
  
  weightR <- weight(cdfR,parWeight)
  
  X <- R*(weightR*Wu + (1-weightR)/(1/Wl))
  Y <- R*(weightR*(1-Wu) + (1-weightR)/(1/(1-Wl)))
  
  
  
  if(onlyXY){
    res <- cbind(X,Y)
    
    if(PLOT){
      par(mfrow=c(1,1))
      plot(X,Y,pch=20,asp=1,xlab="X",ylab="Y",main="Simulated (X,Y) sample") 
    }
  } else{
    XY <- cbind(X,Y)
    res <- list(XY=XY,R=R,Wl=Wl,Wu=Wu,weightR=weightR)
    
    if(PLOT){
      par(mfrow=c(3,4))
      ### first row of plots
      hist(R,breaks=50,freq=FALSE, xlab = "R",main="Histogram of R")
      lines(seq(0,max(R),length=1000),dR(seq(0,max(R),length=1000),parR),col="red",main="Histogram of R")
      hist(1/Wl,breaks=50,xlim=c(range(1/Wl)[1],range(1/Wl)[2]),freq=FALSE, xlab="L",main="Histogram of L")
      lines(seq(range(1/Wl)[1],range(1/Wl)[2],length=1000),dL(seq(range(1/Wl)[1],range(1/Wl)[2],length=1000),parWl),col="red",main="Histogram of L")
      hist(Wu,breaks=50,xlim=c(0,1),freq=FALSE, xlab="U",main="Histogram of U")
      lines(seq(0,1,length=1000),dWu(seq(0,1,length=1000),parWu),col="red",main="Histogram of U")
      plot(seq(0,1,length=1000),weight(seq(0,1,length=1000),parWeight),type="l",xlab="F(R)",ylab="w(F(R))",main="Weight function")
      
      ### second row of plots
      marg.thr <- 0.5
      plot(X,Y,pch=20,asp=1,xlab="X1",ylab="X2",main="Simulated X sample",cex.main=1)
      points(X[X+Y>quantile(X+Y,0.98) & rank(X)/(n+1) > marg.thr & rank(Y)/(n+1) > marg.thr],Y[X+Y>quantile(X+Y,0.98) & rank(X)/(n+1) > marg.thr & rank(Y)/(n+1) > marg.thr],pch=20,cex=1.1,col="red")
      points(X[(1/(X+Y))>quantile(1/(X+Y),0.98) & rank(1/X)/(n+1) > marg.thr & rank(1/Y)/(n+1) > marg.thr],Y[(1/(X+Y))>quantile((1/(X+Y)),0.98) & rank(1/X)/(n+1) > marg.thr & rank(1/Y)/(n+1) > marg.thr],pch=20,cex=1.1,col="blue")
      abline(h=0,v=0,col="lightgrey")
      plot(rank(X)/(n+1),rank(Y)/(n+1),pch=20,asp=1,xlab="X1",ylab="X2",main="Sim X sample on Unif(0,1) scale",cex.main=.95)
      points(rank(X)[(1/(X+Y))>quantile((1/(X+Y)),0.98) & rank(1/X)/(n+1) > marg.thr & rank(1/Y)/(n+1) > marg.thr]/(n+1),rank(Y)[(1/(X+Y))>quantile(1/(X+Y),0.98) & rank(1/X)/(n+1) > marg.thr & rank(1/Y)/(n+1) > marg.thr]/(n+1),pch=20,cex=1.1,col="blue")
      points(rank(X)[X+Y>quantile(X+Y,0.98) & rank(X)/(n+1) > marg.thr & rank(Y)/(n+1) > marg.thr]/(n+1),rank(Y)[X+Y>quantile(X+Y,0.98) & rank(X)/(n+1) > marg.thr & rank(Y)/(n+1) > marg.thr]/(n+1),pch=20,cex=1.1,col="red")
      abline(h=c(0,1),v=c(0,1),col="lightgrey")
      hist((1/(X+Y)),breaks=50,freq=FALSE,xlab="1/|X|", main="Histogram of 1/|X|",cex.main=1)
      abline(v=quantile((1/(X+Y)),0.98),col="blue")
      hist(X+Y,breaks=50,freq=FALSE,xlab="|X|",main="Histogram of |X|",cex.main=1)
      lines(seq(0,max(R),length=1000),dR(seq(0,max(R),length=1000),parR),col="red")
      
      ### third row of plots
      hist(((1/X)/(1/(X+Y)))[(1/(X+Y))>quantile(1/(X+Y),0.98) & rank(1/X)/(n+1) > marg.thr & rank(1/Y)/(n+1) > marg.thr],breaks=50,xlim=c(range(1/Wl)[1],range(1/Wl)[2]),freq=FALSE,xlab="(1/X)/(1/|X|) given 1/|X|>u",main="Hist of (1/X)/(1/|X|) given 1/|X|>u",cex.main=.98, cex.lab=1)
      lines(seq(range(1/Wl)[1],range(1/Wl)[2],length=1000),dL(seq(range(1/Wl)[1],range(1/Wl)[2],length=1000),parWl),col="red",main="Histogram of L")
      hist((X/(X+Y))[X+Y>quantile(X+Y,0.98) & rank(X)/(n+1) > marg.thr & rank(Y)/(n+1) > marg.thr],breaks=50,xlim=c(0,1),freq=FALSE,xlab= "X/|X| given |X|>u", main="Hist of X/|X| given |X|>u",cex.main=1, cex.lab=1)
      lines(seq(0,1,length=1000),dWu(seq(0,1,length=1000),parWu),col="red",main="Histogram of U")
    }
  }
  
  return( res )
}
