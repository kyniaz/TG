#############################
# Incorporando covariáveis; #
#############################

source('importacao_dados.R')

######### 

dados_tg |> head()

dados_tg$idade_cat = ifelse(dados_tg$age_at_diagnosis < 14600, 
                                         "young", ifelse(dados_tg$age_at_diagnosis < 25000,
                                          "not so young", "old"))

y1 = ifelse(dados_tg$idade_cat == "young", 1, 0)
y2 = ifelse(dados_tg$idade_cat == "not so young", 1, 0)
y3 = ifelse(dados_tg$trat == 'Info ausente', 1, 0)
y4 = ifelse(dados_tg$trat == 'Medicamento', 1, 0)
y5 = ifelse(dados_tg$trat == 'Medicamento e Radioterapia', 1, 0)
y6 = ifelse(dados_tg$trat == 'Radioterapia', 1, 0)

t = dados_tg$tempo/365
d = dados_tg$cens

log_veros = function(par){
  b_par =  exp(par[2] + par[3]*y1 + par[4]*y2)
  eta = par[5] + par[6]*y1 + par[7]*y2
  mu = 1/(1 + exp(-eta))
  
  l = 0
  for(i in 1:nrow(dados_filtrados_final)){
    li = d[i]*log(1 - mu[i]) + d[i]*dgompertz(t[i], a = par[1], b = b_par[i], log_opt = T) +
      (1- d[i])*log(mu[i]  + (1-mu[i])*(1 - pgompertz(t[i], a = par[1], b = b_par[i])))
    l = l + li
  }
  return(l)
}

param = c(1, 1, 0.05, 0.05, 0.3, 0, 0)

teste = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(0.0001, -Inf, -Inf,-Inf, -Inf, -Inf, -Inf, -Inf),
              upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf))

teste$par

plot(survfit(Surv(tempo/365, cens) ~ idade_cat, data = dados_tg), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

curve(1/(1 + exp(-teste$par[5])) + (1-1/(1 + exp(-teste$par[5])))*(1-pgompertz(x, teste$par[1], exp(teste$par[2]))), col = 'red', add = T)

curve(1/(1 + exp(-(teste$par[5] + teste$par[6]))) + (1-1/(1 + exp(-(teste$par[5] + teste$par[6]))))*(1-pgompertz(x, teste$par[1], exp(teste$par[2] + teste$par[3]))) , col = 'blue', add = T)

curve(1/(1 + exp(-(teste$par[5] + teste$par[7]))) + (1-1/(1 + exp(-(teste$par[5] + teste$par[7]))))*(1-pgompertz(x, teste$par[1], exp(teste$par[2] + teste$par[4]))) , col = 'orange', add = T)

plot(survfit(Surv(tempo/365, cens) ~ trat, data = dados_tg), 
     xlab = "Anos", 
     ylab = "Sobrevivência")
