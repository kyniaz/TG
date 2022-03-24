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
      out[i] = log(a) - log(b) + x[i]/b -a*(expm1(x[i]/b));
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
NumericVector dgompertz_vec(NumericVector x, double a, NumericVector b, bool log_opt = false) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (log_opt == true){
    for(int i = 0; i < n; ++i) {
      out[i] = log(a/b[i]) + log(x[i]/b[i]) -a*(expm1(x[i]/b[i]));
    }
  }
  else {
    for(int i = 0; i < n; ++i) {
      out[i] = (a/b[i])*exp(x[i]/b[i])*exp(-a*(expm1(x[i]/b[i])));
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericVector pgompertz(NumericVector x, double a, double b, bool log_opt = false, bool lower_tail = true) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else if (lower_tail == true){
    if (log_opt == true){
      for(int i = 0; i < n; ++i) {
        out[i] = log(1 - exp(-a*expm1(x[i]/b)));
      }
    }
    else {
      for(int i = 0; i < n; ++i) {
        out[i] = 1 - exp(-a*(expm1(x[i]/b)));
      }
    }
  }
  else {
    if (log_opt == false){
      for(int i = 0; i < n; ++i) {
        out[i] = -a*expm1(x[i]/b);
      }
    }
    else {
      for(int i = 0; i < n; ++i) {
        out[i] = exp(-a*(expm1(x[i]/b)));
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pgompertz_vec(NumericVector x, double a, NumericVector b, bool log_opt = false) {
  int n = x.size();
  NumericVector out(n);
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (log_opt == true){
    for(int i = 0; i < n; ++i) {
      out[i] = log(1 - exp(-a*expm1(x[i]/b[i])));
    }
  }
  else {
    for(int i = 0; i < n; ++i) {
      out[i] = 1 - exp(-a*(expm1(x[i]/b[i])));
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
        out[i] =  b*log1p(-(1/a)*log1p(-p[i]));
    }
  }
  return out;
}

// [[Rcpp::export]]
double metropolis_gibbs_1(double a, double b, NumericVector t, NumericVector cens, int BB, double start){
  NumericVector chain(BB);
  int n = t.size();
  long double veros_original;
  long double veros_proposta;
  NumericVector dist_original(n);
  NumericVector dist_proposta(n);
  NumericVector sobr_original(n);
  NumericVector sobr_proposta(n);  
  long double prop;
  long double lratio;
  
  chain[0] = start;

  for(int k = 1; k < BB; ++k){
    veros_original = 0;
    veros_proposta = 0;
    
    prop = R::rgamma(chain[k-1], 2) + 0.000001;

    dist_original = dgompertz(t, a, b, true);
    dist_proposta = dgompertz(t, prop, b, true);
    sobr_original = sgompertz(t, a, b, true);
    sobr_proposta = sgompertz(t, prop, b, true);
    
    for(int l = 0; l < n; ++l){
      veros_original = veros_original + cens[l]*dist_original[l] + (1-cens[l])*sobr_original[l];
      veros_proposta = veros_proposta + cens[l]*dist_proposta[l] + (1-cens[l])*sobr_proposta[l];
    }
    
    veros_original = veros_original + R::dgamma(chain[k-1], 1, 2,  true) + R::dgamma( b, 4, 1, true);
    veros_proposta = veros_proposta + R::dgamma(prop, 1, 2, true) + R::dgamma( b, 4, 1, true);
    
    lratio = veros_proposta - veros_original + R::dgamma(prop, chain[k-1], 2, true) - R::dgamma(chain[k-1], prop, 2, true); 
    
    if(log(R::runif(0.0, 1.0)) <= lratio){
      chain[k] = prop;
    }
    else{
      chain[k] = chain[k-1];
    }
  }
  return(chain[chain.size() - 1]);
}

// [[Rcpp::export]]
double metropolis_gibbs_2(double a, double b, NumericVector t, NumericVector cens, int BB, double start){
  NumericVector chain(BB);
  int n = t.size();
  long double prop;
  long double veros_original;
  long double veros_proposta;
  NumericVector dist_original(n);
  NumericVector dist_proposta(n);
  NumericVector sobr_original(n);
  NumericVector sobr_proposta(n);  
  long double lratio;

  chain[0] = start;
  
  for(int k = 1; k < BB; ++k){
    veros_original = 0;
    veros_proposta = 0;
    
    prop = R::rgamma(chain[k-1], 1) + 0.000001;
  
    dist_original = dgompertz(t, a, b, true);
    dist_proposta = dgompertz(t, a, prop, true);
    sobr_original = sgompertz(t, a, b, true);
    sobr_proposta = sgompertz(t, a, prop, true);
    
    for(int l = 0; l < n; ++l){
        veros_original = veros_original + cens[l]*dist_original[l] + (1-cens[l])*sobr_original[l];
        veros_proposta = veros_proposta + cens[l]*dist_proposta[l] + (1-cens[l])*sobr_proposta[l];
    }
    
    veros_original = veros_original + R::dgamma(a, 1, 1,  true) + R::dgamma(b, 4, 1, true);
    veros_proposta = veros_proposta + R::dgamma(a, 1, 1, true) + R::dgamma(prop, 4, 1, true);
    
    lratio = veros_proposta - veros_original + R::dgamma(prop, chain[k-1], 1, true) - R::dgamma(chain[k-1], prop, 1, true); 
    
    if(lratio > log(R::runif(0.0, 1.0))){
        chain[k] = prop;
      }
      else{
        chain[k] = chain[k-1];
      }
  }
  return(chain[chain.size() - 1]);
}

// [[Rcpp::export]]
NumericVector metropolis_gibbs_C(NumericVector t, NumericVector cens, int B, int BB, double start_a, double start_b) {
  NumericMatrix out(2, B); //Matriz de parâmetros
  int dim;
  
  out(0, 0) = start_a;
  out(1, 0) = start_b;
  
  for(int i = 1; i < B; ++i) {
    dim = trunc(R::runif(0.0, 2));
    
    switch (dim)
    {
    case 0:
      out(0, i) = metropolis_gibbs_1(out(0, i - 1), out(1, i - 1), t, cens, BB, start_a);
      out(1, i) = out(1, i - 1);
      break;
    case 1:
      out(0, i) = out(0, i - 1);
      out(1, i) = metropolis_gibbs_2(out(0, i - 1), out(1, i - 1), t, cens, BB, start_b);
      break;
    }
  }
  return(out);
}