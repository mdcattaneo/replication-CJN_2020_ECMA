#############################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Apr-24-2017
##  This version: June-3-2020
##  Function File
#############################################
library(Rcpp)
sourceCpp("main_maxscore_cpp.cpp")

## setting seed
SeedFunction <- function() {
    args = commandArgs(trailingOnly=TRUE)
    
    # error handling
    if (length(args)==0) { args <- list("1") }
    
    # the argument is
    cat(paste("The argument is", args[[1]], "\n", sep=" "))
    
    # generate seed sequence
    set.seed(42)
    Seeds <- sample(1:50000, 50000, replace=FALSE)
    
    # set the seed
    cat(paste("The current seed is", toString(Seeds[as.integer(args[[1]])]), "\n", sep=" "))
    set.seed(Seeds[as.integer(args[[1]])])
    
}



## generating observations
dgp <- function(N,beta,MODEL){
    X1 <- rnorm(N, mean = 0, sd = 1)
    X2 <- rnorm(N, mean = 1, sd = 1)
  
    if (MODEL == 1) U <- rlogis(N,location=0, scale= sqrt(3)/pi)/sqrt(2)
    if (MODEL == 2) U <- rt(N, df=3)/sqrt(3)
    if (MODEL == 3){
        Z <- X1 + X2
        V <- rlogis(N, location = 0, scale = sqrt(3)/pi )
        U <- 0.25*(1 + 2*Z^2 + Z^4)*V
    }
  
    Y <- as.numeric( X1 + X2*beta + U >= 0 )
  
    return(list(y = Y, x1 = X1, x2 = X2))
}

## Objective function of max score estimator
Mn <- function(OBS,beta){
    mean( (2*OBS$y - 1)*(OBS$x1 + OBS$x2*beta > 0) )
}

######################################################
## Plug-in estimator for Hessian: \hat{H}_{MS}
######################################################
## kernel function
K1 <- function(x) -x*dnorm(x)

hat.H.ms <- function(OBS, beta, bw){
    mean( (2*OBS$y - 1)*K1( (OBS$x1 + OBS$x2*beta)/bw )/bw^2*OBS$x2^2 )
}

#############################################################################
## Rule-of-thumb bandwidth for plug-in estimator and Numerical derivative
## Using heteroscedastic probit
#############################################################################
ll <- function(y,x1,x2,par){
    z <- x1 + x2*par[1]
    v <- par[2] + par[3]*z + par[4]*z^2 + par[5]*z^3 + par[6]*z^4
    t1 <- y*log(pnorm(z/v) )
    t2 <- (1-y)*log(1-pnorm(z/v))
    t1[is.nan(t1)] <- -Inf
    t2[is.nan(t2)] <- -Inf
    min( -mean( t1 + t2 ), 1e3)
}

ROT <- function(DATA){
    y <- DATA$y
    x1 <- DATA$x1
    x2 <- DATA$x2
    
    n <- length(y)
    tmp <- optim( c(1,1,rep(0,4)), ll, y = y, x1 = x1, x2 = x2)
    para <- tmp$par
    v <- para[2]
    v1 <- para[3]
    v2 <- 2*para[4]
    mu1 <- mean(x1)
    sigma1 <- sd(x1)
    
    F1 <- dnorm(0)/v
    b1 <- F1*( (x2*para[1] + mu1)^2/sigma1^2 - 1)/sigma1^3*dnorm( (x2*para[1] + mu1)/sigma1 )
    F3 <- dnorm(0)*(v2/v^2 - 2*v1^2/v^3) - dnorm(0)/v^3
    b2 <- 1/3*F3*dnorm( (x2*para[1] + mu1)/sigma1 )/sigma1
    
    B.ker <- mean( x2^2*b1 + x2^2*b2)*(-3)
    B.nd <- mean( x2^4*b1 + x2^4*b2 )/2
    
    V.ker <- mean( x2^2*dnorm( (x2*para[1] + mu1)/sigma1 )/sigma1 )*0.1410474
    V.nd <- mean( abs(x2)*dnorm( (x2*para[1] + mu1)/sigma1 )/sigma1 )*2
    
    ker.h <- (3*V.ker/4/B.ker^2)^(1/7)*n^(-1/7)
    nd.h <- (3*V.nd/4/B.nd^2)^(1/7)*n^(-1/7)
    
    return( list( bw.kernel = ker.h, bw.nd = nd.h) )
    
}

###################################################################
## Numericative derivative estimator for Hessian: \hat{H}
###################################################################
hat.H <- function(OBS, beta, eps){
    ( Mn(OBS, beta + eps) + Mn(OBS, beta - eps) - 2*Mn(OBS, beta) )/eps^2
}



#########################################################################################
## MERGE OUTPUT FILES
#########################################################################################
merge.output = function(cpus=100,models=1:3){
  
    for (m in models) {
        files = NULL;
        for (s in 1:cpus) {
            filename = paste0("output/parts/N_ms_m",m,"_",s,".txt")
            if (file.exists(filename)) files = c(files, filename)
        }
        if ( length(filename) > 0 ){
            write.table(do.call("rbind", lapply(files, function(x) read.table(x))), file = paste0("output/N_ms_m",m,".txt") )
            message("\nModel ",m," done -- ",length(files)," Files read.")
        } else message("\nModel ",m," not available.")
    }
}
