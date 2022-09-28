library('survival')
library('ggplot2')
library('flexsurv')
####### Sem censura ----

survfit_gompertz = function(y, d, dados = NULL) {
  log_like = function(par) {
    a = par[1]
    b = par[2]
    return(sum(d*dgompertz(y, a, b, T) + (1-d)*(pgompertz(y, a, b, T, lower.tail = F))))
  }
  
  start = c(1,1)
  
  param = optim(start, log_like, control = list(fnscale = -1, maxit = 500),
                method="L-BFGS-B", lower = c(0.0001, 0.0001), upper = c(Inf,Inf))
  
  print(param$value)
  return(param$par)
}


set.seed(154)

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5)
)

lista_tempos = list()
for(i in 1:5){
  tempos = rgompertz(n[i], 0.5, 2)
  
  lista_tempos[[i]] = tempos
  
  mod = survfit_gompertz(tempos, rep(1, n[i]))
  estimativas[i,c(1,2,3)] = c(n[i], 
                              mod[1], 
                              mod[2])
}

estimativas

xtable::xtable(estimativas, digits = 3) |> print(include.rownames = F)

dados_simulados = data.frame(tempos, rep(1, max(n)))
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

dados1 = plot_dados
print(xtable::xtable(plot_dados |> tail(n = 10), digits = 3), include.rownames = F)
print(xtable::xtable(plot_dados |> head(n = 10), digits = 3), include.rownames = F)

######Para cada conjunto, plotar a curva encima do km;

xplot = seq(0.01, 2.5, 0.01)

camadas = list()

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a[i], estimativas$b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

freq_g = ggplot() +
  geom_step(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = plot_dados, size = 1) +
  labs(x = 'Tempo') +
  theme_minimal() + 
  geom_line(aes(x=x, y=d, colour = "b", linetype = "b"), data = camadas[[1]], size = 1) + 
  geom_line(aes(x=x, y=d, colour = "c", linetype = "c"), data = camadas[[2]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "d", linetype = "d"), data = camadas[[3]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "e", linetype = "e"), data = camadas[[4]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "f", linetype = "f"), data = camadas[[5]], size = 1) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(#"deeppink2",
                       #"darkorchid3",
                       "dodgerblue",
                       #"gray12",
                       "springgreen3",
                       "lightsteelblue",
                       "tomato",
                       "gold",
                       "royalblue4"),
                     labels = c("Kaplan-Meier", 
                                "n = 25", 
                                "n = 50",
                                "n = 100",
                                "n = 200",
                                "n = 500")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("Kaplan-Meier", 
                                   "n = 25", 
                                   "n = 50",
                                   "n = 100",
                                   "n = 200",
                                   "n = 500")) 

ggsave(filename = 'figuras/curvas_estimadas_sem_cens.pdf', units = 'in',
       width = 7, height = 5)
# Com censura ----

x_teste = seq(0.001, 2, 0.001)

plot(x_teste, dgompertz(x_teste, 1, 5), type = "l", col = "blue")
lines(x_teste, dgompertz(x_teste, 1, 1), col = "red")

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5),
  cens = numeric(5)
)

set.seed(154)

lista_tempos = list()
for(i in 1:5){
  t1 = rgompertz(n[i], 1, 5)
  t2 = rgompertz(n[i], 1, 1)
  
  eventos = t1 < t2
  
  tempos = c(t1[eventos], t2[eventos == F])
  
  lista_tempos[[i]] = tempos
  
  mod = survfit_gompertz(tempos, eventos)
  estimativas[i,c(1,2,3)] = c(n[i], 
                              mod[1], 
                              mod[2])
  
  estimativas[i, 4] = 1 - sum(eventos)/n[i]
}

estimativas

xtable::xtable(estimativas, digits = 3) |> print(include.rownames = F)

dados_simulados = data.frame(tempos, eventos)
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

dados2 = plot_dados
print(xtable::xtable(plot_dados |> tail(n = 10), digits = 3), include.rownames = F)
print(xtable::xtable(plot_dados |> head(n = 10), digits = 3), include.rownames = F)

######Para cada conjunto, plotar a curva encima do km;

xplot = seq(0.01, 1, 0.01)
camadas = list()

par(mfrow = c(1,1))

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a[i], 
                                                estimativas$b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

freq_cc = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = plot_dados, size = 1) +
  labs(x = 'Tempo') +
  theme_minimal() + 
  geom_line(aes(x=x, y=d, colour = "b", linetype = "b"), data = camadas[[1]], size = 1) + 
  geom_line(aes(x=x, y=d, colour = "c", linetype = "c"), data = camadas[[2]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "d", linetype = "d"), data = camadas[[3]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "e", linetype = "e"), data = camadas[[4]], size = 1) +
  geom_line(aes(x=x, y=d, colour = "f", linetype = "f"), data = camadas[[5]], size = 1) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(#"deeppink2",
                       #"darkorchid3",
                       "dodgerblue",
                       #"gray12",
                       "springgreen3",
                       "lightsteelblue",
                       "tomato",
                       "gold",
                       "royalblue4"),
                     labels = c("Kaplan-Meier", 
                                "n = 25", 
                                "n = 50",
                                "n = 100",
                                "n = 200",
                                "n = 500")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("Kaplan-Meier", 
                                   "n = 25", 
                                   "n = 50",
                                   "n = 100",
                                   "n = 200",
                                   "n = 500")) 
freq_cc

ggsave(filename = 'figuras/curvas_estimadas_com_cens.pdf', units = 'in',
       width = 7, height = 5)

fim = cbind(dados1, dados2)

print(xtable::xtable(fim |> tail(n = 10), digits = 3), include.rownames = F)
print(xtable::xtable(fim |> head(n = 10), digits = 3), include.rownames = F)
