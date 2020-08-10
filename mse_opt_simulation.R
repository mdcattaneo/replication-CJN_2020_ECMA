##################################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Apr-6-2017
##  This version: June-3-2020
##  MSE-optimal tuning parameters via simulations
###################################################
rm(list=ls())

####################################
## Functions
####################################
library(Rcpp)

sourceCpp("Hhat.cpp")
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



###############################################
## Simulation setup
###############################################

H0 <- unlist( read.table("H0.txt") )

models <- 1:3
beta0 <- 1

bw.grid <- seq(0.01, 2 ,0.01)
bwlen <- length(bw.grid)

eps.grid <- seq(0.01, 2, 0.01)
epslen <- length(eps.grid)

n <- 1000
S <- 5e4

seed.start <- 99
seed.inc <- 6

MSE.ms <- matrix(NA, nrow = S, ncol = bwlen)
MSE.nd <- matrix(NA, nrow = S, ncol = epslen)

opt.tuning.para <- matrix(NA, nrow = 2, ncol = length(models))
#### Monte Carlo
time <- Sys.time()

for (m in models){
    for (s in 1:S){
        seed <- seed.start + seed.inc*(s-1)
        set.seed(seed)
        
        obs <- dgp(n,beta0,m)
        
        MSE.ms[s,] <- ( H_ms(obs$y, obs$x1, obs$x2, bw.grid, beta0) - H0[m] )^2
        MSE.nd[s,] <- ( H(obs$y, obs$x1, obs$x2, eps.grid, beta0) - H0[m] )^2
        
        
        if (s%%100==0){
            cat(paste("\nSimulations Completed:",s,
                      "Diff:", round(difftime(Sys.time(),time,units="mins"),2),"mins"));
            time <- Sys.time();
        }
    }
    ####################################
    ## Find optimal tuning parameters
    ####################################
    ave.MSE.ms <- colMeans(MSE.ms)
    ave.MSE.nd <- colMeans(MSE.nd)
    
    opt.tuning.para[1,m] <- bw.grid[ which.min(ave.MSE.ms) ]
    opt.tuning.para[2,m] <- eps.grid[ which.min(ave.MSE.nd) ]
}


write.table(opt.tuning.para, file = "mseopt_bw.txt")
