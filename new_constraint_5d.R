library(fields)
library(mvtnorm)
library(flexclust)
library(geoR)
library(Matrix)
library(truncnorm)
library(MCMCpack)
library(splines)
library(tmg)

b <- read.csv("C:/Users/why94/Dropbox/materialcodes/b.csv")
A <- read.csv("C:/Users/why94/Dropbox/materialcodes/A.csv")
A$X <- NULL
b$X <- NULL
b <- as.matrix(b)
A <- as.matrix(A)
A <- Matrix(A)

n <- 3000

sigma.sq <- 5
tau.sq1 <- 0.5
phi1 <- 1
tau.sq2 <- 0.5
phi2 <- 3
tau.sq3 <- 0.5
phi3 <- 5

D <- as.matrix(dist(b))
R1 <- exp(-phi1*D)
R2 <- exp(-phi2*D)
R3 <- exp(-phi3*D)
cc <- kronecker(R,diag(3))

w1 <- rnorm(n, 0, sigma.sq*R1)
w2 <- rnorm(n, 0, sigma.sq*R2)
w3 <- rnorm(n, 0, sigma.sq*R3)
y1 <- rnorm(n, w1, sqrt(tau.sq))
y2 <- rnorm(n, w2, sqrt(tau.sq))
y3 <- rnorm(n, w3, sqrt(tau.sq))
y <- as.vector(t(cbind(y1,y2,y3)))

res <- A %*% y
# ind <- sample(1:l,25)
# A[ind,n] <- A[ind,n] - b[ind] / y[n] 
A[,3*n] <- A[,3*n] - res / y[3*n] 
yy <- rep(0,n)

predbbasis <- bs(b[,1],df=10)
for (j in c(2,3,5,6)){
  predbbasis1 <- bs(b[,j],df=10)
  predbbasis <- cbind(predbbasis,predbbasis1)
}
predbbasis <- cbind(1,predbbasis)

predbasis <- matrix(0,3*n,3*ncol(predbbasis))
for (i in 1:n){
  predbasis[3*i-2,1:ncol(predbbasis)] <- predbbasis[i,]
  predbasis[3*i-1,1:ncol(predbbasis)+ncol(predbbasis)] <- predbbasis[i,]
  predbasis[3*i,1:ncol(predbbasis)+2*ncol(predbbasis)] <- predbbasis[i,]
}


Bayes_Kriging<-function(y,s,basis,A,
                        sp=NULL,Xp=NULL,
                        df=4,b_var=diag(dim(A)[1]),
                        sd_beta=1000,can_df = 4,
                        iters=1000,burn=100){
  
  
  # Bookkeeping
  
  p        <- ncol(basis)
  
  
  predictions<-!is.null(sp) 
  np <- 1
  if(predictions){
    np <- nrow(sp)*3
    predbbasis <- bs(sp[,1],df=10)
    for (j in c(2,3,5,6)){
      predbbasis1 <- bs(sp[,j],df=10)
      predbbasis <- cbind(predbbasis,predbbasis1)
    }
    predbbasis <- cbind(1,predbbasis)
    
    predbasis <- matrix(0,3*np,3*ncol(predbbasis))
    for (i in 1:n){
      predbasis[3*i-2,1:ncol(predbbasis)] <- predbbasis[i,]
      predbasis[3*i-1,1:ncol(predbbasis)+ncol(predbbasis)] <- predbbasis[i,]
      predbasis[3*i,1:ncol(predbbasis)+2*ncol(predbbasis)] <- predbbasis[i,]
    }
  }
  
  
  # Initial values
  # bbb <- matrix(0,n,2*n)
  # for (i in 1:n){
  #   bbb[i,(2*i-1):(2*i)] <- s[i,]
  # }
  X <- A%*%basis
  beta    <- rep(1,dim(basis)[2])
  fff <- beta
  sigma2 <- cov(matrix(y-basis%*%fff,n,3))
  Sigma <- A%*%kronecker(diag(n),sigma2)%*%t(A)
  
  Siginv  <- solve(Sigma + diag(dim(Sigma)[1]))
  
  
  # Keep track of stuff
  
  keep.beta <- matrix(0,iters,p)
  keeper.sigma2 <-   list()
  pred  <- pred_norm    <- matrix(0,iters,np)
  
  
  #colnames(keepers)   <- c("sigma2","range","gamma2","epsilon")
  colnames(keep.beta) <- colnames(X)
  
  # GO!!!
  
  for(i in 1:iters){
    if(i %% 1 ==0){print(i)}
    ##############################################:
    #####       beta(Gibbs)    #######:
    ##############################################:
    tXS  <- t(X)%*% Siginv
    VVV1  <- tXS%*%X + diag(p)/sd_beta^2
    VVV <- solve(VVV1)
    MMM  <- VVV%*%tXS%*%y
    Q <- list()
    Q[[1]] <- list(t(basis)%*%basis, rep(0,p), -n)
    beta <- as.vector(rtmg(1,as.matrix(VVV1),as.vector(VVV1%*%MMM),fff,f = basis,
                           g = rep(0,dim(basis)[1]),q=Q))
    fff <- beta
    
    
    #############################################:
    ####          sigma2 (Gibbs)        #######:
    #############################################:
    
    R      <- y- X%*%beta
    a      <- n + df
    b      <- R%*%t(R) + b_var
    Sigma <- riwish(a,b_var)
    Siginv  <- solve(Sigma + diag(dim(Sigma)[1]))
    
    
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.beta[i,]  <- beta
    keeper.sigma2[[i]]  <- Sigma
    
    
    
    ##############################################:
    #####           PREDICTIONS            #######:
    ##############################################:
    
    if(i>burn & predictions){
      pred[i,] <- predbasis%*%beta
      pred_norm[i,] <- pred[i,] / sqrt(sum(pred[i,]^2))
    }
    
  }   
  list(beta=keep.beta,sigma2=keeper.sigma2,pred=pred, pred_norm=pred_norm)}


try <- Bayes_Kriging(yy,b,predbasis,A ,iters = 1000)
