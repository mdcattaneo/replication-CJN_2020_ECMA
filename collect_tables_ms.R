#############################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Apr-24-2017
##  This version: May-10-202
##  Constructing tables
#############################################
rm(list=ls())
library(Hmisc)

models <- 1:3

beta0 <- 1

out <- matrix(NA, nrow = 10, ncol = 3*3) #ncol is average bandwidth, coverage probability, ave CI length x 3 DGPs

for (m in models){
    results <- read.table( paste0("output/N_ms_m",m,".txt") )
    
    ## lower and upper ends of CI
    ci.low <- results[,"bhat"] - results[,"q.975"]
    ci.upp <- results[,"bhat"] - results[,"q.025"]
    ## Whether CI contains the true parameter value
    test.tmp <- ( ci.low <= beta0) & ( beta0 <= ci.upp )
    ## Length of CI
    ci.tmp <- ci.upp - ci.low
    
    ## average tuning parameters
    out[,3*m-2] <- by(results[,"h.e"], results[,"label"], mean)
    ## frequency with which CI contains the true parameter
    out[,3*m-1] <- by(test.tmp, results[,"label"], mean)
    ## average length of CI
    out[,3*m] <- by(ci.tmp, results[,"label"], mean)
}    

    
output <- latex(round(out,3), file = "table_maxscore.txt",
                   append=FALSE, table.env=FALSE, center="none", title="",
                   n.cgroup=c(3, 3, 3),
                   cgroup=c("DGP 1", "DGP 2", "DGP 3"),
                   colheads= rep( c("$h,\\epsilon$","Coverage", "Length"),3 ),
                   n.rgroup = c(1,3,3,3),
                   rgroup = c("Standard","m-out-of-n","Plug-in: $\\tilde{\\mathbf{H}}^{\\mathtt{MS}}_n$", 
                              "Num Deriv: $\\tilde{\\mathbf{H}}^{\\mathtt{ND}}_n$"),
                   rowname=c("", "$m = \\lceil n^{1/2} \\rceil$", "$m = \\lceil n^{2/3} \\rceil$",
                             "$m = \\lceil n^{4/5} \\rceil$", 
                             "$h_{ \\mathtt{MSE} } $","$h_{\\mathtt{AMSE} }$", "$\\hat{h}_{\\mathtt{AMSE} }$", 
                             "$\\epsilon_{ \\mathtt{MSE} } $", "$\\epsilon_{\\mathtt{AMSE} }$", "$\\hat{\\epsilon}_{\\mathtt{AMSE} }$") )
    


