library('survival')
library('ggplot2')

####### Sem censura ----


set.seed(154)

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5)
)

#Frequentista

#amostra

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

xtable::xtable(estimativas, digits = 3)

dados_simulados = data.frame(tempos, rep(1, max(n)))
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

print(xtable(plot_dados |> tail(n = 20), digits = 3), include.rownames = F)

estimativas

######Para cada conjunto, plotar a curva encima do km;
xplot = seq(0.01, 5, 0.01)

#fo
camadas = list()

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a[i], estimativas$b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

freq_g = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a"), data = plot_dados, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "c"), data = camadas[[2]]) +
  geom_line(aes(x=x, y=d, colour = "d"), data = camadas[[3]]) +
  geom_line(aes(x=x, y=d, colour = "e"), data = camadas[[4]]) +
  geom_line(aes(x=x, y=d, colour = "f"), data = camadas[[5]]) +
  scale_color_manual(name = "",
                     values = c(
                       "brown2",
                       "dodgerblue",
                       "springgreen3",
                       "#CBEB55",
                       "#4DCDF5",
                       "gold"),
                     labels = c("Kaplan-Meier","n = 25",
                                "n = 50",
                                "n = 100",
                                "n = 200",
                                "n = 500")) +
  theme(legend.position = 'bottom')
# Com censura ----

x_teste = seq(0.001, 6, 0.001)

plot(x_teste, dgompertz(x_teste, 1, 1.35), type = "l", col = "blue")
lines(x_teste, dgompertz(x_teste, 1, 2), col = "red")

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5),
  cens = numeric(5)
)

set.seed(154)
#Frequentista

#amostra

lista_tempos = list()
for(i in 1:5){
  t1 = rgompertz(n[i], 1, 1.35)
  t2 = rgompertz(n[i], 0.3, 2)
  
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

dados_simulados = data.frame(tempos, eventos)
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

#print(xtable(plot_dados |> tail(n = 20), digits = 3), include.rownames = F)

estimativas

######Para cada conjunto, plotar a curva encima do km;

xplot = seq(0.01, 3, 0.01)
camadas = list()

par(mfrow = c(1,1))

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a[i], 
                                                estimativas$b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

freq_cc = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a"), data = plot_dados, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "c"), data = camadas[[2]]) +
  geom_line(aes(x=x, y=d, colour = "d"), data = camadas[[3]]) +
  geom_line(aes(x=x, y=d, colour = "e"), data = camadas[[4]]) +
  geom_line(aes(x=x, y=d, colour = "f"), data = camadas[[5]]) +
  scale_color_manual(name = "",
                     values = c(
                       "brown2",
                       "dodgerblue",
                       "springgreen3",
                       "#CBEB55",
                       "#4DCDF5",
                       "gold"),
                     labels = c("Kaplan-Meier","n = 25",
                                "n = 50",
                                "n = 100",
                                "n = 200",
                                "n = 500")) +
  theme(legend.position = 'bottom')

library(xtable)

print(xtable(estimativas, digits = 3),include.rownames = F)
