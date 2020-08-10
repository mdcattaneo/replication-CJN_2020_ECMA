#################################################
## Cattaneo, Jansson and Nagasawa (2020, ECMA)
## Replication Files
#################################################

R files

-DGP_constants.R
Compute H_0 (Hessian matrix, here it is a scalar as we focus on one-dimensional settings)

-mse_opt_simulation.R
Using H_0 computed by DGP_constants.R, conduct Monte Carlos to find MSE-optimal tuning parameters to estimate H_0

-AMSE.R
Compute Asymptotic MSE and the optimal tuning parameters based on the leading term of the AMSE.

-main_maxscore.R
This is the main R file that conducts simulations for bootstrap-based inference.
The code was run in batch mode i.e., using R CMD BATCH

-main_function_maxscore.R
This file is called by main_maxscore.R and defines functions used within main_function.R

-collect_tables_ms.R
This file compiles the output from main_maxscore.R and creates the table in the paper. 

C++ files
-Hhat.cpp
Called within mse_opt_simulation.R to compute the estimator for H_0

-main_maxscore_cpp.cpp
This file is called by main_function_maxscore.R and computes the maximum score estimator and other objects.

Text files
-H0.txt
Contains the values of H_0

-mseopt_bw.txt
Contains the MSE optimal tuning parameters computed using simulations

-AMSEopt_bw.txt
Contains the optimal tuning parameters based on the leading term of the AMSE.