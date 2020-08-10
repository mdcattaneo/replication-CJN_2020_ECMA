#########################################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Apr-6-2017
##  This version: May-10-2020
##  Computing Approximate MSE optimal tuning parameters
##########################################################
rm(list=ls())
n <- 1000
models <- 1:3
AMSE.tuning.para <- matrix( NA, nrow = 2, ncol = length(models))

##########################################
## Plug-in Estimator
##########################################
fz0 <- function(x) exp(-x^2/2)/sqrt(2*pi)
fz1 <- function(x) -x*exp(-x^2/2)/sqrt(2*pi)
fz2 <- function(x) (x^2-1)*exp(-x^2/2)/sqrt(2*pi)

Ku3 <- function(x) (-x*exp(-x^2/2)/sqrt(2*pi))*x^3
Ksqr <- function(x) (-x*exp(-x^2/2)/sqrt(2*pi))^2

k3 <- -3 #integrate(Ku3,-Inf, Inf)$value
kappa <- integrate(Ksqr,-Inf,Inf)$value



v.intgrd <- function(x) x^4*fz0(x)*dnorm(x,mean=1,sd=1)

V <- integrate(v.intgrd, -Inf, Inf)$value*kappa

##################################
##    DGP1
###################################
s <- sqrt(3)/pi/sqrt(2)

Fe1 <- 1/4/s
Fe2 <- 0
Fe3 <- -1/8/s^3


B1.intgrd <- function(x) -Fe1*fz2(x)*x^2*dnorm(x,1,1)
B3.intgrd <- function(x) -Fe3*fz0(x)*x^2*dnorm(x,1,1)/3

B1 <- integrate(B1.intgrd, -Inf,Inf)$value*k3
B3 <- integrate(B3.intgrd, -Inf,Inf)$value*k3

Bsqr <- (B1 + B3)^2

AMSE.tuning.para[1,1] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)

##################################
##    DGP2
###################################
C <- 2/pi

Fe1 <- C
Fe2 <- 0
Fe3 <- -4*C

b1.intgrd <- function(x) -Fe1*fz2(x)*x^2*dnorm(x,1,1)
b3.intgrd <- function(x) -Fe3*fz0(x)*x^2*dnorm(x,1,1)/3

B1 <- integrate(b1.intgrd, -Inf,Inf)$value*k3
B3 <- integrate(b3.intgrd, -Inf,Inf)$value*k3

Bsqr <- (B1+B3)^2

AMSE.tuning.para[1,2] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)

##################################
##    DGP3
###################################
s <- sqrt(3)/pi
Fe1 <- 4*1/s/4
Fe2 <- 0
Fe3 <- 64*(-1/8/s^3) - 32*(1/s/4)

b1.intgrd <- function(x) -Fe1*fz2(x)*x^2*dnorm(x,1,1)
b3.intgrd <- function(x) -Fe3*fz0(x)*x^2*dnorm(x,1,1)/3

B1 <- integrate(b1.intgrd,-Inf,Inf)$value*k3
B3 <- integrate(b3.intgrd,-Inf,Inf)$value*k3

Bsqr <- (B1 + B3)^2

AMSE.tuning.para[1,3] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)

####################################
## Numerical Derivative
####################################
v.intgrd <- function(x) abs(x)*fz0(x)*dnorm(x,mean=1,sd=1)
V <- integrate(v.intgrd, -Inf, Inf)$value*2


##################################
##    DGP1
###################################
s <- sqrt(3)/pi/sqrt(2)

Fe1 <- 1/4/s
Fe2 <- 0
Fe3 <- -1/8/s^3


B1.intgrd <- function(x) fz2(x)*x^4*dnorm(x,1,1)
B3.intgrd <- function(x) fz0(x)*x^4*dnorm(x,1,1)

B1 <- -integrate(B1.intgrd, -Inf,Inf)$value*(-1/2)
B3 <- -integrate(B3.intgrd, -Inf,Inf)$value*(-1/6)

Bsqr <- (B1*Fe1 + B3*Fe3)^2

AMSE.tuning.para[2,1] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)

##################################
##    DGP2
###################################
C <- 2/pi

Fe1 <- C
Fe2 <- 0
Fe3 <- -4*C

Bsqr <- (B1*Fe1 + B3*Fe3)^2

AMSE.tuning.para[2,2] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)
##################################
##    DGP3
###################################
s <- sqrt(3)/pi
Fe1 <- 4*1/s/4
Fe2 <- 0
Fe3 <- 64*(-1/8/s^3) - 32*(1/s/4)

Bsqr <- (B1*Fe1 + B3*Fe3)^2

AMSE.tuning.para[2,3] <- (3*V/4/Bsqr)^(1/7)*n^(-1/7)


####################
write.table(AMSE.tuning.para, file = "AMSEopt_bw.txt")