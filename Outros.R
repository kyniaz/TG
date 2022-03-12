######Para cada conjunto, plotar a curva encima do km;
xplot = seq(0.01, max(t), 0.01)

camadas = list()

nova_camada = data.frame(xplot, teste$par[3] + (1-teste$par[3])*(1-pgompertz(xplot, teste$par[1], teste$par[2])))
colnames(nova_camada) = c('x','d')
camadas[[1]] = nova_camada

nova_camada = data.frame(xplot, teste2$par[3]^(pgompertz(xplot, a = teste2$par[1], b = teste2$par[2])))
colnames(nova_camada) = c('x','d')
camadas[[2]] = nova_camada

emv_bc = plot_sim + geom_line(aes(x=x, y=d, colour = "a"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[2]]) +
  scale_color_manual(name = "Estimativa",
                     values = c(
                       "springgreen4",
                       "brown2"),
                     labels = c("Com Mistura", "Sem Mitura")) +
  theme(legend.position = 'bottom')

ggsave('figuras_TG/EMV_bc.pdf', emv_bc, units = 'in', width = 7, height = 5)

plot(density(a))
median(a)
median(b)
median(prob)

#curve(median(prob) + (1-median(prob))*(1-pgompertz(x, median(a), median(b))), col = 'red', add = T)
curve(median(prob)^(pgompertz(x, a = median(a), b = median(b))), col = 'blue', add = T)

#Gerar da gompertz;
rgompertz()

hazard = function(x, a, b){
  return(dgompertz(x,a,b)/(1-pgompertz(x,a,b)))
}

library('ggplot2')

lista_a_bs = list()

lista_a_bs[[1]] = c(-3,1)
lista_a_bs[[2]] = c(-1.5,0.5)
lista_a_bs[[3]] = c(0.25,0.5)
#lista_a_bs[[3]] = c(1,0.5)
lista_a_bs[[4]] = c(0.5,0.5)
lista_a_bs[[5]] = c(1,0.25)
# lista_a_bs[[6]] = c(1,0.75)
# lista_a_bs[[7]] = c(1,1.25)
lista_a_bs[[6]] = c(0,2)

x = seq(0, 3, 0.01)

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
                     labels = c("a = -3; b = 1", 
                                "a = -1.5; b = 0.5", 
                                #"1,0.5",
                                "a = 0.25; b = 0.5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0.5; b = 0.5",
                                "a = 1; b = 0.25",
                                "a = 0; b = 2"))

########## densidade

x = seq(0, 3, 0.01)

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
                     labels = c("a = -3; b = 1", 
                                "a = -1.5; b = 0.5", 
                                #"1,0.5",
                                "a = 0.25; b = 0.5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0.5; b = 0.5",
                                "a = 1; b = 0.25",
                                "a = 0; b = 2"))

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
                     labels = c("a = -3; b = 1", 
                                "a = -1.5; b = 0.5", 
                                #"1,0.5",
                                "a = 0.25; b = 0.5",
                                #"1,0.25",
                                #"1,0.75",
                                "a = 0.5; b = 0.5",
                                "a = 1; b = 0.25",
                                "a = 0; b = 2"))

library('gridExtra')

plos = grid.arrange(haz, dens, ncol = 2)

ggsave('haz_dens.pdf', plos, units = c('in'), height = 5, width = 10)

#Simula??o
#write.csv(dados_simulados, file = "simulado.csv")

n_curados = rgompertz(100, 1, 2)
curados = rgompertz(500, 2, 1)

eventos = c(numeric(100) + 1, numeric(500))

tempos = c(n_curados, curados)
# 
dados_simulados = data.frame(tempos, eventos)

kaplan_meier_s = survfit(Surv(tempos, eventos) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

plot_sim = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'brown2') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(x = Tempo, y = Sobrevivência), data = plot_dados[plot_dados$Evento == 1,], shape = 3, color = 'royalblue', size = 0.5) +
  labs(x = 'Tempo') +
  ylim(0, 1) + 
  theme_minimal() 


n = c(10, 50, 100, 200, 600)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5)
)

for(i in 1:5){
  simulado_sample = dados_simulados[sample(1:600, n[i]),]
  mod = flexsurvreg(Surv(tempos, eventos)~ 1, data = simulado_sample, dist = 'gompertz')
  estimativas[i,] = c(n[i],  mod$coefficients[1], exp(mod$coefficients)[2])
}

estimativas

######Para cada conjunto, plotar a curva encima do km;
xplot = seq(0.01, 1.6, 0.01)

#fo
camadas = list()

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a[i], estimativas$b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

plot_sim + geom_line(aes(x=x, y=d, colour = "a"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[2]]) +
  geom_line(aes(x=x, y=d, colour = "c"), data = camadas[[3]]) +
  geom_line(aes(x=x, y=d, colour = "d"), data = camadas[[4]]) +
  geom_line(aes(x=x, y=d, colour = "e"), data = camadas[[5]]) +
  scale_color_manual(name = "Tamanho Amostral",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "gold",
                       "orange"),
                     labels = c("10","50","100","200","600")) +
  theme(legend.position = 'bottom')