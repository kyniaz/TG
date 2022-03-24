library('Rcpp')
sourceCpp(file = 'gompertz.cpp')
library('survival')
library('ggplot2')

survfit_gompertz = function(y, d, dados = NULL) {
  log_like = function(par) {
    a = par[1]
    b = par[2]
    return(sum(d*dgompertz(y, a, b, log_opt = T) + (1-d)*(sgompertz(y, a, b, log_opt = T))))
  }
  
  start = c(1,1)
  
  param = optim(start, log_like, control = list(fnscale = -1, maxit = 500),
                method="L-BFGS-B", lower = c(0.0001, 0.0001), upper = c(Inf,Inf))
  
  return(param$par)
}


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
                     labels = c("Kaplan-Meier","Tamanho Amostral: 25",
                                "Tamostra Amostra 50",
                                "Tamanho Amostral 100",
                                "Tamanho Amostral 200",
                                "Tamanho Amostral 500")) +
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
                     labels = c("Kaplan-Meier","Tamanho Amostral = 25",
                                "Tamanho Amostral = 50",
                                "Tamanho Amostral = 100",
                                "Tamanho Amostral = 200",
                                "Tamanho Amostral = 500")) +
  theme(legend.position = 'bottom')

library(xtable)
print(xtable(estimativas, digits = 3),include.rownames = F)

# Risco, Dens e Sobr ----

hazard = function(x, a, b){
  #return(dgompertz(x,a,b, log_opt = T) - (sgompertz(x,a,b, log_opt = T)))
  (a/b)*exp(x/b)
}

lista_a_bs = list()

lista_a_bs[[1]] = c(3,0.5)
lista_a_bs[[2]] = c(2,0.5)
lista_a_bs[[3]] = c(1,0.5)
#lista_a_bs[[3]] = c(1,0.5)
lista_a_bs[[4]] = c(0.5,0.25)
lista_a_bs[[5]] = c(1, 2)
# lista_a_bs[[6]] = c(1,0.75)
# lista_a_bs[[7]] = c(1,1.25)
lista_a_bs[[6]] = c(2,2)


x = seq(0.001, 1, 0.001)

figuras = lapply(1:6, function(n) hazard(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)


haz = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6"),
            size = 1)+
  # geom_line(data=dados_risco, 
  #           mapping=aes(x=x, y=C7, colour = "7"),
  #           size = 1)+
  # geom_line(data=dados_risco, 
  #           mapping=aes(x=x, y=C8, colour = "8"),
  #           size = 1)+
  # geom_line(data=dados_risco, 
  #           mapping=aes(x=x, y=C9, colour = "9"),
  #           size = 1)+
  theme_minimal()+
  labs(y = 'h(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
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
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                #"1,0.5",
                                "a = 1; b = 0,5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2"))

ggsave('figuras/haz_gompertz.pdf', haz, units = 'in', width = 7, height = 5)
########## densidade

figuras = lapply(1:6, function(n) dgompertz(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)


dens = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6"),
            size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C7, colour = "7"),
  #          size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C8, colour = "8"),
  #          size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C9, colour = "9"),
  #          size = 1)+
  theme_minimal()+
  labs(y = 'f(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
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
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                #"1,0.5",
                                "a = 1; b = 0,5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2"))

ggsave('figuras_TG/dens_gompertz.pdf', dens, units = 'in', width = 7, height = 5)

#Sobrevivência

x = seq(0, 3, 0.01)

figuras = lapply(1:6, function(n) 1 - pgompertz(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)


survi = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6"),
            size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C7, colour = "7"),
  #          size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C8, colour = "8"),
  #          size = 1)+
  #geom_line(data=dados_risco, 
  #          mapping=aes(x=x, y=C9, colour = "9"),
  #          size = 1)+
  theme_minimal()+
  labs(y = 'S(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
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
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                #"1,0.5",
                                "a = 1; b = 0,5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2"))

ggsave('figuras_TG/surv_gompertz.pdf', survi, units = 'in', width = 7, height = 5)
