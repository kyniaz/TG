source('importacao_dados.R')

######################
#TTT plot com censura#
######################

tempo = dados_filtrados_final$tempo
censura = dados_filtrados_final$cens

o=order(tempo)
t=tempo[o]
cens=censura[o]


n=length(t)
r=sum(cens)

j=1
TF=numeric()
MON=numeric()
I=numeric()
TTT=numeric()
Fi=numeric()
F_var=numeric()
S=numeric()

TF[j]=0
MON[j]=0
F_var[j]=0
S[j]=1
TTT[j]=0
i=1

while(i<(n+1)){
  if(cens[i]==1){
    j=j+1
    TF[j]=t[i]
    NI=n-i+1
    I=((n+1)-MON[j-1])/(1+NI)
    MON[j]=MON[j-1]+I
    F_var[j]=MON[j]/n
    S[j]=1-F_var[j]
  }
  i=i+1
}

TF[r+2]=t[n]
F_var[r+2]=1
TTT[1]=0

for(j in 2:(r+2)){
  TTT[j]=TTT[j-1]+n*S[j-1]*(TF[j]-TF[j-1])
}

for(j in 1:(r+2)){
  Fi[j]=TTT[j]/TTT[r+2]
}

ggplot(data = data.frame(F_var, Fi), aes(x = F_var, y = Fi))+
  geom_point() +
  theme_minimal() + 
  labs(x = "F(t)", y = "TTT com censuras") +
  geom_abline(slope=1, intercept=0) +
  geom_line() + 
  xlim(0,1) + ylim(0,1)


kaplan_meier = survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final)

plot_dados = data.frame(kaplan_meier$time, kaplan_meier$surv, kaplan_meier$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'steelblue') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #geom_point(aes( x = Tempo, y = Sobrevivência), data = plot_dados[plot_dados$Evento == 1,], shape = 3, color = 'brown2') +
  labs(x = 'Tempo em anos') +
  ylim(0, 1) + 
  theme_minimal() 

par(mfrow = c(2,2))

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

min_optim = 10^(-5)
############### GOMPERTZ ----

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dgompertz(t, a = par[1], b = par[2], log_opt = T) + (1-d)*log(par[3] + (1-par[3])*(1-pgompertz(t, a = par[1], b = par[2])))))
}

param = c(1,1,0.5)

mix_model = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

mix_model$par

curve(mix_model$par[3] + (1-mix_model$par[3])*(1-pgompertz(x, mix_model$par[1], mix_model$par[2])), col = 'red', add = T)


#Veros ponderada
-log(mix_model$par[3])*mix_model$value

log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dgompertz(t, a = par[1], b = par[2], log_opt = T) + log(par[3])*pgompertz(t, a = par[1], b = par[2]) ))
}

param2 = c(1,1,0.3)

no_mix_model = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

no_mix_model$par 


curve((no_mix_model$par[3])^(pgompertz(x, a = no_mix_model$par[1], b = no_mix_model$par[2])), col = 'blue', add = T)

############### GAMMA -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dgamma(t, par[1], par[2], log = T) + 
               (1-d)*log(par[3] + (1-par[3])*(1-pgamma(t, par[1], par[2])))))
}

param = c(1,1,0.5)


mix_model = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

mix_model$par

curve(mix_model$par[3] + (1-mix_model$par[3])*(1-pgamma(x, mix_model$par[1], mix_model$par[2])), col = 'red', add = T)


-log(mix_model$par[3])*mix_model$value

log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dgamma(t, par[1], par[2], log = T) + 
               log(par[3])*pgamma(t, par[1], par[2]) ))
}

param2 = c(1,1,0.3)

no_mix_model = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

no_mix_model$par 


curve((no_mix_model$par[3])^(pgamma(x, no_mix_model$par[1], no_mix_model$par[2])), col = 'blue', add = T)

########## Exponencial -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[2])) + d*dexp(t, par[1], log = T) + 
               (1-d)*log(par[2] + (1-par[2])*(1-pexp(t, par[1])))))
}

param = c(0.01,0.5)

mix_model = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

mix_model$par

curve(mix_model$par[2] + (1-mix_model$par[2])*(1-pexp(x, mix_model$par[1])), col = 'red', add = T)

log_veros2 = function(par){
  return(sum(d*log(-log(par[2])) + d*dexp(t, par[1], log = T) + 
               log(par[2])*pexp(t,par[1]) ))
}

param2 = c(0.01,0.3)

no_mix_model = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

no_mix_model$par 

curve((no_mix_model$par[2])^(pexp(x, no_mix_model$par[1])), col = 'blue', add = T)

########## WEIBULL -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dweibull(t, par[1], par[2], log = T) + 
               (1-d)*log(par[3] + (1-par[3])*(1-pweibull(t, par[1], par[2])))))
}

param = c(1,1,0.5)


mix_model = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

mix_model$par

curve(mix_model$par[3] + (1-mix_model$par[3])*(1-pweibull(x, mix_model$par[1], mix_model$par[2])), col = 'red', add = T)


log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dweibull(t, par[1], par[2], log = T) + 
               log(par[3])*pweibull(t, par[1], par[2]) ))
}

param2 = c(1,1,0.3)

no_mix_model = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(min_optim,min_optim,min_optim), upper = c(Inf,Inf,0.9999))

no_mix_model$par 


curve((no_mix_model$par[3])^(pweibull(x, no_mix_model$par[1], no_mix_model$par[2])), col = 'blue', add = T)
