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


double hat_V_ms(NumericVector y, NumericVector x1, NumericVector x2, double bw, double beta){
    int n = y.size();
    double out = 0.0;
    for (int i = 0; i < n; i++){
        out += (2.0*y[i] - 1.0)*ker( (x1[i] + x2[i]*beta)/bw )/pow(bw, 2)*pow(x2[i], 2);
    }
    out = out/n;
    return(out);
}


//[[Rcpp::export]]
NumericVector V_ms(NumericVector y, NumericVector x1, NumericVector x2, NumericVector bws, double beta){
    int lenh = bws.size();
    NumericVector out(lenh,0.0);
    
    for (int i =0; i < lenh; i++){
        out[i] = hat_V_ms(y, x1, x2, bws[i], beta);
    }
    return(out);
}


//[[Rcpp::export]]
NumericVector V(NumericVector y, NumericVector x1, NumericVector x2, NumericVector grds, double beta){
    int len = grds.size();
    NumericVector out(len, 0.0);
    
    for (int i = 0; i < len; i++){
        out[i] = ( score_fnc(y,x1,x2,beta + grds[i]) + score_fnc(y,x1,x2,beta - grds[i]) 
                   - 2*score_fnc(y,x1,x2,beta) )/pow(grds[i], 2);
    }
    return(out);
}