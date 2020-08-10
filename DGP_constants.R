#############################################
##  Bootstrapping Maximum Score Estimator
##  Cattaneo Jansson Nagasawa
##  Mar-6-2017
##  This version: June-3-2020
##  Compute H_0
#############################################
rm(list=ls())
H <- rep(NA,3)
## DGP: Y = \I\{ x_1 +x_2 + u\geq 0 \}
## (x_1,x_2) =_d Normal( (0,1), I_2 ) where I_2 is the 2x2 identity matrix
## The distribution of u varies across the 3 DGPs

## H_0= 2E[|x_2|^2 f_{u|x_1x_2}(0|-x_2,x_2) f_{x_1|x_2}(-x_2|x_2) ]
## For DGP 1 & 2, u is independent of x's and for DGP 3, u|x_1=-x_2, x_2=_d 0.25*v, where v has the logistic distribution with mean 0 and variance 1
## Then, for the 3 DGPs, H_0 = 2 f_{u}(0) E[|x_2|^2 f_{x_1|x_2}(-x_2|x_2) ]
## Also, x_1 is independent of x_2 so f_{x_1|x_2}=f_{x_1}.
## Let C = E[|x_2|^2 f_{x_1}(-x_2) ]
## and using normality,
C <- -3*exp(-1/4)/8/sqrt(pi) 

##################################
##    DGP 1: logistic
###################################
s <- sqrt(3)/pi
fu <- 1/s/4*sqrt(2)

H[1] <- 2*fu*C
##################################
##    DGP 2: t3 distribution
###################################
fu <-  2/pi

H[2] <- 2*fu*C
##################################
##    DGP 3: heteroscedastic
###################################
s <- sqrt(3)/pi/4
fu <- 1/s/4

H[3] <- 2*fu*C

write.table(H, file="H0.txt")