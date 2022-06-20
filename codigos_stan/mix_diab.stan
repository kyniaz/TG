functions {
  real dgompertz_lpdf(real x, real a, real b) {
      //return log(b) + a*x - (b/a)*(expm1(a*x));
      return log(a) - log(b) + x/b -a*(expm1(x/b));
  }
  real sgompertz_lpdf(real x, real a, real b) {
      //return log(exp(-(b/a)*expm1(a*x)));
      return log(exp(-a*expm1(x/b)));
  }
  real log_veros_mix_lpdf(real x, int d, real a, real b, real theta) {
      //return d*(log(1 - theta)) + d*dgompertz_lpdf(x|a,b) + (1-d)*log(theta + (1 - theta)*exp(-(b/a)*(expm1(a*x))));
      return d*(log(1 - theta)) + d*dgompertz_lpdf(x|a,b) + (1-d)*log(theta + (1 - theta)*exp(-a*(expm1(x/b))));
  }
}

data {
  int<lower=0> N;  
  real<lower=0> T[N];             
  int D[N]; 
}

parameters {
  real<lower=0.000001> a; // parametros do modelo.
  real<lower=0.000001> b;
  real<lower=0,upper=1> theta;
}

model {
  a ~ gamma(40, 1);        // prioris
  b ~ gamma(50, 0.5); 
  theta ~ beta(1, 1);
  for (k in 1:N) {
    T[k] ~ log_veros_mix_lpdf(D[k], a, b, theta);
  }
}
