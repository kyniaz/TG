library(rstan)

dados_tg = read.csv('dados_tg.csv', header = T)

acoisa = "
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
  real log_veros_nomix_lpdf(real x, int d, real a, real b, real theta) {
      return d*(log(-log(theta))) + d*dgompertz_lpdf(x|a,b) + log(theta)*(1 - exp(-a*(expm1(x/b))));
  }
}

data {
  int<lower = 0> N;  
  real<lower=0> T[N];               // number of trials
  int D[N];   // success on trial n
}

parameters {
  real<lower=0.0001> a; // chance of success
  real<lower=0.0001> b;
  real<lower=0,upper=1> theta;
}

model {
  a ~ gamma(1, 1);        // prior
  b ~ gamma(6, 1); 
  theta ~ beta(1, 1);
  for (k in 1:N) {
    T[k] ~ log_veros_nomix_lpdf(D[k], a, b, theta);
  }
}
"

fit = stan(model_code=acoisa,
           data=list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
           iter=5000, chains=2, cores=2)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

acoisa2 = "
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
  real log_veros_nomix_lpdf(real x, int d, real a, real b, real theta) {
      return d*(log(-log(theta))) + d*dgompertz_lpdf(x|a,b) + log(theta)*(1 - exp(-a*(expm1(x/b))));
  }
  real log_veros_mix_lpdf(real x, int d, real a, real b, real theta) {
      return d*(log(1 - theta)) + d*dgompertz_lpdf(x|a,b) + (1-d)*log(theta + (1 - theta)*exp(-a*(expm1(x/b))));
  }
}

data {
  int<lower = 0> N;  
  real<lower=0> T[N];               // number of trials
  int D[N];   // success on trial n
}

parameters {
  real<lower=0.0001> a; // chance of success
  real<lower=0.0001> b;
  real<lower=0,upper=1> theta;
}

model {
  a ~ gamma(1, 1);        // prior
  b ~ gamma(6, 1); 
  theta ~ beta(1, 1);
  for (k in 1:N) {
    T[k] ~ log_veros_mix_lpdf(D[k], a, b, theta);
  }
}
"

fit2 = stan(model_code=acoisa2,
           data=list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
           iter=5000, chains=2, cores=2)

print(fit2, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit2)


fit_summary2 = summary(fit2)

fit_summary2$summary

fit_summary2$summary[,1]

x = dados_tg$tempo/365

kaplan_meier_s = survfit(Surv(tempo/365, cens) ~ 1, data = dados_tg)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "Kaplan-Meier"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y= fit_summary2$summary[,1]['theta'] + (1 - fit_summary2$summary[,1]['theta'])*pgompertz(x, fit_summary2$summary[,1]['a'], 
                                            fit_summary2$summary[,1]['b'], lower.tail = F),
                colour = "Com Mix"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= fit_summary$summary[,1]['theta']^pgompertz(x, fit_summary$summary[,1]['a'], 
                                                                   fit_summary$summary[,1]['b']),
                colour = "No Mix"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom')
