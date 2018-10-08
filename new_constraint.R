library(fields)
library(mvtnorm)
library(flexclust)
library(geoR)

library(truncnorm)
library(MCMCpack)
library(splines)
library(tmg)

fxx <- function(s){
  y <- rep(0,length(s))
  for (i in 1:length(s)){
    if(s[i]<2*pi){y[i] <- (sin(s[i]))}
    else if(s[i] < 4*pi){
      y[i] <- (-5*sin(s[i] / 2))}
    else{y[i] <- ((s[i] - 4*pi)^2 / 10)}
  }
  y <- ifelse(s*y >= 1,y,1 / s)
  #y <- sqrt(length(s)) * y / sqrt(sum(y^2))
  return(y)}


n <- 100
l <- 30
x <- seq(0.5,20,length.out = n)
y <- fxx(x)
plot(fxx(x),type = "l")

sp <- x


set.seed(0531)
A <- matrix(0,l,n)
# for (i in 1:20){
#   A[i,(3*i-2):(3*i)] <- rnorm(3,0,1)
# }
# for (i in 21:l){
#   A[i,((i-20)*4-3):((i-20)*4) + 60] <- rnorm(4,0,1)
# }
for (i in 1:l){
  ind <- sample(1:n,10)
  A[i,ind] <- rnorm(10)
}


b <- A %*% y
# ind <- sample(1:l,25)
# A[ind,n] <- A[ind,n] - b[ind] / y[n] 
A[,n] <- A[,n] - b / y[n] 

bbasis <- bs(x,df=10)
bbasis <- cbind(1,bbasis)


Bayes_Krige<-function(y,s,X,
                      sp=NULL,
                      a_var=.01,b_var=.01,
                      sd_beta=1000,
                      iters=1000,burn=100){
  
  
  # Bookkeeping
  
  n        <- length(y)
  p        <- ncol(X)
  
  
  predictions<-!is.null(sp) 
  np <- 1
  if(predictions){
    predbasis <- bs(sp,df = 10)
    predbasis <- cbind(1,predbasis)
    np <- length(sp)
  }
  
  # Initial values
  X <- A %*% X
  #beta <- c(0,0,0)
  beta    <- rep(1,p)
  sigma2  <- var(y - X%*%beta)[1,1]
  Sigma   <- A %*% t(A) * sigma2
  Siginv  <- solve(Sigma + diag(n))
  
  # Keep track of stuff
  
  keep.beta <- matrix(0,iters,p)
  keepers   <- matrix(0,iters,2)
  pred  <- pred_norm    <- matrix(0,iters,np)
  
  colnames(keepers)   <- c("sigma2","range")
  colnames(keep.beta) <- colnames(X)
  
  # GO!!!
  fff <- beta
  for(i in 1:iters){
    
    ##############################################:
    #####       MEAN PARAMETERS (Gibbs)    #######:
    ##############################################:
    if (i %% 10 == 0) {print(i)}
    tXS  <- t(X)%*%Siginv
    VVV1 <- tXS%*%X + diag(p)/sd_beta^2
    VVV  <- solve(VVV1)
    MMM  <- VVV%*%tXS%*%y
    # fff <- lm((apply(cbind(1/s,st-1),1,max)+1.2)~bbasis-1)$coefficients
    # beta <- t(rtmg(1,VVV1,as.vector(VVV1%*%MMM),fff,f=rbind(bbasis,-bbasis),
    #              g=c(-apply(cbind(1/s,st-1),1,max),st+2.01)))
    Q <- list()
    Q[[1]] <- list(t(bbasis)%*%bbasis, rep(0,p), -n)
    beta <- as.vector(rtmg(1,VVV1,as.vector(VVV1%*%MMM),fff,f = bbasis,
                           g = rep(0,dim(bbasis)[1]),q=Q))
    fff <- beta
    ##############################################:
    #####          VARIANCE (Gibbs)        #######:
    ##############################################:
    
    R      <- y-X%*%beta
    a      <- n/2+a_var
    b      <- t(R)%*%Siginv%*%R/2+b_var
    sigma2 <- 1/rgamma(1,a,b) 
    Sigma   <- A %*% t(A) * sigma2
    Siginv  <- solve(Sigma + diag(n))
    ##############################################:
    #### CORRELATION PARAMETERS (Metropolis) #####:
    ##############################################:
    
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.beta[i,]  <- beta
    keepers[i,1]   <- sigma2
    #keepers[i,2]   <- exp(theta[1])
    
    
    ##############################################:
    #####           PREDICTIONS            #######:
    ##############################################:
    
    if(i>burn & predictions){
      pred[i,] <- predbasis%*%beta
      pred_norm[i,] <- pred[i,] / sqrt(sum(pred[i,]^2))
    }
    
  }   
  
  
  list(beta=keep.beta,keepers=keepers,pred=pred,pred_norm = pred_norm)}

xp <- x
try <- Bayes_Krige(rep(0,30),x,bbasis,xp)



lb <- apply(try$pred[100:1000,],2,quantile,probs = c(0.025,0.975))
plot(x,y,type = "l")
lines(x,colMeans(try$pred_norm[100:1000,]),col = "red")
lines(x,lb[1,],col = "red",lty = 2)
lines(x,lb[2,],col = "red", lty = 2)

ratio <- max(colMeans(try$pred_norm[100:1000,]))/max(y)
y_pred <- colMeans(try$pred_norm[100:1000,]) / ratio

plot(x,y,type = "l")
lines(x,y_pred , type="l",col="red")
lb1 <- apply(try$pred_norm[100:1000,],2,quantile,probs = c(0.025,0.975))
lines(x,lb1[1,] / ratio,col = "red",lty = 2)
lines(x,lb1[2,] / ratio,col = "red", lty = 2)
