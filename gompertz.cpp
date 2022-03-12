#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dgompertz(NumericVector x, double a, double b, bool log_opt = false) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else if (log_opt == true){
    for(int i = 0; i < n; ++i) {
      out[i] = log(a/b) + log(x[i]/b) -a*(expm1(x[i]/b));
    }
  }
  else {
    for(int i = 0; i < n; ++i) {
      out[i] = (a/b)*exp(x[i]/b)*exp(-a*(expm1(x[i]/b)));
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericVector pgompertz(NumericVector x, double a, double b, bool log_opt = false) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else if (log_opt == true){
    for(int i = 0; i < n; ++i) {
      out[i] = log(1 - exp(-a*expm1(x[i]/b)));
    }
  }
  else {
    for(int i = 0; i < n; ++i) {
      out[i] = 1 - exp(-a*(expm1(x[i]/b)));
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector sgompertz(NumericVector x, double a, double b, bool log_opt = false) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else if (log_opt == true){
    for(int i = 0; i < n; ++i) {
      out[i] = -a*expm1(x[i]/b);
    }
  }
  else {
    for(int i = 0; i < n; ++i) {
      out[i] = exp(-a*(expm1(x[i]/b)));
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector qgompertz(NumericVector x, double a, double b) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else {
    for(int i = 0; i < n; ++i) {
      if ((x[i] <= 1) & ( x[i] > 0)){
        out[i] = b*log(1 - (1/a)*log(1-x[i]));
      }
      else {
        out[i] = R_NaN;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector rgompertz(double n, double a, double b) {
  NumericVector p = Rcpp::runif(n, 0.0, 1.0);
  NumericVector out(n);
  
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else {
    for(int i = 0; i < n; ++i) {
        out[i] =  b*log(1 - (1/a)*log(p[i]));
    }
  }
  return out;
}