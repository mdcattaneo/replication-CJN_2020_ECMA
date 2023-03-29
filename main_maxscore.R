#############################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Apr-24-2017
##  This version: Mar-28-2020
##  Main File
#############################################
source("main_function_maxscore.R")

## SEED SETUP FOR PARALLELIZATION
SeedFunction()

ags <- commandArgs(trailingOnly = TRUE)
# error handling
if (length(ags)==0) { ags <- list("1") }
ag <- ags[[1]]
#########################################################################################
## MONTECARLO SETUP
#########################################################################################

ncpus <- 20
## sample size
n <- 1000
## Number of bootstrap replications
B <- 2000
## Total number of Monet Carlo iterations
S <- 2000%/%ncpus



## population \beta
beta0 <- 1
models <- 1:3
## grids over which the objective function is maximized
grids <- c( seq(-1,0,0.01), seq(0.002,2,0.002), seq(2.01,3,0.01) )

## subsample size for m-out-of-n bootstrap
subsample <- c( ceiling(n^(1/2) ), ceiling(n^(2/3) ), ceiling(n^(4/5) )  )

## loading MSE optimal bandwidth calculated via simulation
tuning.para <- read.table("mseopt_bw.txt")
## Simulated MSE-optimal bandwidth for \hat{H}^{\mathtt{MS}}: to distinuish from the numerical derivative estimator, we refre to \hat{H}^{\mathtt{MS}} as plug-in estimator
## In the previous version, the \hat{H} was computed over a grid of bandwidths 
# plugin.bw <- rbind(seq(0.7,1.3,0.1)*tuning.para[1,1], seq(0.7,1.3,0.1)*tuning.para[1,2], seq(0.7,1.3,0.1)*tuning.para[1,3])
plugin.bw <- as.numeric(tuning.para[1,])
bw.len <- 1 #previously, it was: bw.len <- length(seq(0.7,1.3,0.1))

## Infeasible AMSE optimal bandwidth
amse.para <- read.table("AMSEopt_bw.txt")
plugin.AMSE.bw <- as.numeric( amse.para[1,] )



## Simulated MSE-optimal tuning parameter for numerical derivative
## In the previous version, the \hat{H} was computed over a grid of bandwidths 
# num.deriv.eps <- rbind(seq(0.7,1.3,0.1)*tuning.para[2,1], seq(0.7,1.3,0.1)*tuning.para[2,2], seq(0.7,1.3,0.1)*tuning.para[2,3])
num.deriv.eps <- as.numeric( tuning.para[2,] )
eps.len <- 1 # previously, it was: eps.len <-length(seq(0.7,1.3,0.1))

## Infeasible AMSE optimal tuning parameter
num.deriv.AMSE <- as.numeric( amse.para[2,] )


num.rows.table <- 10 
# 1 (standard nonparametric bootstrap) + 3 (m-out-of-n) + 2*3 (Simulated-MSE,AMSE, ROT)*(plug-in + numerical derivative) 



col.names <- c("q.025","q.975","h.e", "bhat","label") 
#############################
## Monte Carlo Simulations ##
#############################
time <- Sys.time()

for (m in models){
    cat("\n##########################################")
    cat("\n###         Starting Model ",m,"         ###")
    cat("\n##########################################")
    ## matrix to store results
    out <- matrix(NA, nrow = S*num.rows.table, ncol = 5, dimnames = list(NULL, col.names) )
    ## initialization
    Row <- 1 
    
    for (s in 1:S){
        obs <- dgp(n,beta0,m)
        
        
        ## Maximum score estimator computed from the original sample
        bhat <- maxscore(obs$y, obs$x1, obs$x2, grids)
        
        ## Computing \hat{H}^{\mathtt{MS}} i.e., plug-in estimator
        # H.ms <- rep(NA,bw.len)
        # for (h in 1:bw.len){
        #     H.ms[h] <- hat.H.ms(obs, bhat, plugin.bw[m,h])
        # }
        H.ms <- hat.H.ms(obs, bhat, plugin.bw[m])
        ## AMSE optimal bandwidth (infeasible)
        H.ms.AMSE <- hat.H.ms(obs, bhat, plugin.AMSE.bw[m])
        # ROT bandwidth
        ROTBW <- ROT(obs)
        plugin.ROT.bw <- ROTBW$bw.kernel
        H.ms.ROT <- hat.H.ms(obs, bhat, plugin.ROT.bw)
        
        
        ## computing numerical derivative estimator \hat{H} 
        # H <- rep(NA, eps.len)
        # for (h in 1:eps.len){
        #     H[h] <- hat.H(obs, bhat, num.deriv.eps[m,h])
        # }
        H <- hat.H(obs, bhat, num.deriv.eps[m])
        ## AMSE optimal epsilon (infeasible)
        H.AMSE <- hat.H(obs, bhat, num.deriv.AMSE[m])
        ## ROT epsilon
        num.deriv.ROT.eps <- ROTBW$bw.nd
        H.ROT <- hat.H(obs, bhat, num.deriv.ROT.eps)
        
        
        ##########################
        ## Bootstrap Iterations ##
        ##########################
        ## vectors/matrix to store bootstrap values
        std.boot <- rep(NA, B)
        mn.boot <- matrix(NA, nrow = B, ncol = length(subsample) )
        plugin.MSE <- rep(NA, B)            #plugin.MSE <- matrix(NA, nrow = B, ncol = bw.len)
        plugin.AMSE <- rep(NA, B)
        plugin.ROT <- rep(NA,B)
        nd.MSE <- rep(NA,B)                 #nd.MSE <- matrix(NA, nrow = B, ncol = eps.len)
        nd.AMSE <- rep(NA, B)
        nd.ROT <- rep(NA, B)
          
        for (b in 1:B){
            boot.idx <- sample(n, replace = T)
            boot.y <- obs$y[boot.idx]
            boot.x1 <- obs$x1[boot.idx]
            boot.x2 <- obs$x2[boot.idx]
            boot.obs <- list(y = boot.y, x1 = boot.x1, x2 = boot.x2)
            
            ## standard nonparametric bootstrap 
            std.boot[b] <- maxscore(boot.y, boot.x1, boot.x2, grids)
            
            ## m-out-of-n bootstrap
            for (mn in 1:length(subsample) ){
                mn.boot[b,mn] <- maxscore(boot.y[1:subsample[mn]], boot.x1[1:subsample[mn]], boot.x2[1:subsample[mn]], grids)
            }
            
            ## \hat{H}^{\mathtt{MS}} i.e., plug-in estimator 
            # for (h in 1:bw.len){
            #     plugin.MSE[b,h] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids, H.ms[h], bhat)
            # }
            plugin.MSE[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids, H.ms, bhat)
            plugin.AMSE[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H.ms.AMSE, bhat)
            plugin.ROT[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H.ms.ROT, bhat)
            
            ## numerical differentiation
            # for (h in 1:eps.len){
            #     nd.MSE[b,h] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H[h], bhat)
            # }
            nd.MSE[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H, bhat)
            nd.AMSE[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H.AMSE, bhat)
            nd.ROT[b] <- boot_ms(obs$y, obs$x1, obs$x2, boot.obs$y, boot.obs$x1, boot.obs$x2, grids,H.ROT, bhat)
              
        }
        ## end of bootstrap iterations
        
        ##compute quantiles
        
        for (k in 1:num.rows.table){
            if (k == 1){ #standard nonpara bootstrap
                out[Row, 1:2] <- quantile(std.boot - bhat, probs = c(0.025,0.975), type = 1)
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (2 <= k & k <= 4){ #m-out-of-n bootstrap
                out[Row, 1:2] <- quantile(mn.boot[,k -1] - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 5){ #plug-in \hat{H}, MSE optimal bandwidth grids
                out[Row, 1:2] <- quantile(plugin.MSE - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- plugin.bw[m]
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 6){ #plug-in \hat{H}, AMSE optimal bandwidth
                out[Row, 1:2] <- quantile(plugin.AMSE - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- plugin.AMSE.bw[m]
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 7){ #plug-in \hat{H}, ROT bandwidth
                out[Row, 1:2] <- quantile(plugin.ROT - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- plugin.ROT.bw
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 8){ #numerical derivative \hat{H}, MSE optimal grids
                out[Row, 1:2] <- quantile(nd.MSE - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- num.deriv.eps[m]
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 9){ #numerical derivative \hat{H}, AMSE epsilon
                out[Row, 1:2] <- quantile(nd.AMSE - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- num.deriv.AMSE[m]
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            if (k == 10){ #numerical derivative \hat{H}, ROT epsilon
                out[Row, 1:2] <- quantile(nd.ROT - bhat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 3] <- num.deriv.ROT.eps
                out[Row, 4] <- bhat
                out[Row, 5] <- k
            }
            Row = Row + 1
        }
        cat(paste("\nSimulations Completed:",s,
                  "Diff:", round(difftime(Sys.time(),time,units="mins"),2),"mins"));
        time <- Sys.time();
    }
    
    ##################
    ## Save Results ##
    ##################
    write.table(out, file = paste0("output/parts/N_ms_m",m,"_",ag,".txt") )
}

merge.output(ncpus, models)
