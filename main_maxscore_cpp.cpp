#include <Rcpp.h>

using namespace Rcpp;


double ker(double x){
    double anss = -x/sqrt(2*PI)*exp(-pow(x,2)/2);
    return(anss);
}


double score_fnc(NumericVector y, NumericVector x1, NumericVector x2, double beta){
    int n = y.size();
    double obj = 0.0;
    
    for (int i=0; i < n; i++){
        obj += (2.0*y[i]-1.0)*( x1[i]+x2[i]*beta >= 0 ) ;
    }
    obj = obj/n;
    
    return(obj);
}


//[[Rcpp::export]]
double maxscore(NumericVector y, NumericVector x1, NumericVector x2, NumericVector bgrids){
    int lenb = bgrids.size();
    NumericVector objfnc(lenb,0.0);
    double out;
    
    for (int i = 0; i < lenb; i++){
        objfnc[i] = score_fnc(y,x1,x2,bgrids[i]);
    }
    out = bgrids[which_max(objfnc)];
    return(out);
}



//[[Rcpp::export]]
double boot_ms(NumericVector y, NumericVector x1, NumericVector x2, NumericVector bootY, 
               NumericVector bootX1, NumericVector bootX2, NumericVector grids, double H, double beta0){
    int lenb = grids.size();
    NumericVector objfnc(lenb, 0.0);
    double out;
    
    for (int b = 0; b < lenb; b++){
        objfnc[b] = score_fnc(bootY, bootX1, bootX2, grids[b]) - score_fnc(y, x1, x2, grids[b])
                    + H/2*pow(grids[b] - beta0, 2);
    }
    out = grids[ which_max(objfnc) ];
    return(out);
}