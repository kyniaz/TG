functions {
  real dgompertz_lpdf(real x, real a, real b) {
      return log(b) + a*x - (b/a)*(expm1(a*x));
  }
  real sgompertz_lpdf(real x, real a, real b) {
      return (-b/a)*expm1(a*x);
  }
  real log_veros_lpdf(real x, int d, real a, real b) {
      return d*dgompertz_lpdf(x|a,b) + (1-d)*sgompertz_lpdf(x|a,b);
  }
}

data {
  int<lower = 0> N;  
  real<lower=0> T[N];               // number of times
  int D[N];   // censored
}

parameters {
  real a; // a
  real<lower=0.000001> b; // b
}

model {
  a ~ normal(0, 0.5);        // prior
  b ~ gamma(0.2, 2); 
  for (k in 1:N) {
    T[k] ~ log_veros_lpdf(D[k], a, b);
  }
}