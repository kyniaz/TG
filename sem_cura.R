plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

min_optim = 10^(-5)
############### GOMPERTZ ----

log_veros = function(par){
  return(sum(d*dgompertz(t, a = par[1], b = par[2], log_opt = T) + (1-d)*(1-pgompertz(t, a = par[1], b = par[2]))))
}

param = c(1,1)

mix_model = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

mix_model$par

curve((1-pgompertz(x, mix_model$par[1], mix_model$par[2])), col = 'red', add = T)

