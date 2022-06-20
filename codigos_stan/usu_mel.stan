functions {
  real dgompertz_lpdf(real x, real a, real b) {
      return log(a) - log(b) + x/b -a*(expm1(x/b));
  }
  real sgompertz_lpdf(real x, real a, real b) {
      return log(exp(-a*expm1(x/b)));
  }
  real log_veros_lpdf(real x, int d, real a, real b) {
      return d*dgompertz_lpdf(x|a,b) + (1-d)*sgompertz_lpdf(x|a,b);
  }
}

data {
  int<lower=0> N;  
  real<lower=0> T[N];               
  int D[N];   
}

parameters {
  real<lower=0.000001> a; 
  real<lower=0.000001> b;
}

model {
  a ~ gamma(65, 0.5);      
  b ~ gamma(1000, 1); 
  for (k in 1:N) {
    T[k] ~ log_veros_lpdf(D[k], a, b);
  }
}
