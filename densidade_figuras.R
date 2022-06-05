# Figuras de Risco, Densidade e Sobrevivencia ----
source('gompertz.R')
library('ggplot2')

## Risco instantaneo ----
hazard = function(x, a, b){
  (a/b)*exp(x/b)
}

lista_a_bs = vector(mode = "list", length = 6)

lista_a_bs[[1]] = c(3,0.5)
lista_a_bs[[2]] = c(2,0.5)
lista_a_bs[[3]] = c(1,0.5)
lista_a_bs[[4]] = c(0.5,0.25)
lista_a_bs[[5]] = c(1, 2)
lista_a_bs[[6]] = c(2,2)

x = seq(0.001, 1, 0.001)

figuras = lapply(1:6, function(n) hazard(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)


haz = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1", linetype = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2", linetype = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3", linetype = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4", linetype = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5", linetype = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6", linetype = "6"),
            size = 1)+
  theme_minimal()+
  labs(y = 'h(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "tomato",
                       "gold",
                       "royalblue4"),
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                "a = 1; b = 0,5",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("a = 3; b = 0,5", 
                                   "a = 2; b = 0,5", 
                                   "a = 1; b = 0,5",
                                   "a = 0,5; b = 0,25",
                                   "a = 1; b = 1",
                                   "a = 2; b = 2"))

ggsave('figuras/haz_gompertz.pdf', haz, units = 'in', width = 7, height = 5)

## Densidade ----

figuras = lapply(1:6, function(n) dgompertz(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)

dens = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1", linetype = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2", linetype = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3", linetype = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4", linetype = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5", linetype = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6", linetype = "6"),
            size = 1)+
  theme_minimal()+
  labs(y = 'f(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "tomato",
                       "gold",
                       "royalblue4"),
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                "a = 1; b = 0,5",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("a = 3; b = 0,5", 
                                   "a = 2; b = 0,5", 
                                   "a = 1; b = 0,5",
                                   "a = 0,5; b = 0,25",
                                   "a = 1; b = 1",
                                   "a = 2; b = 2"))

ggsave('figuras/dens_gompertz.pdf', dens, units = 'in', width = 7, height = 5)

## Sobrevivência ----

x = seq(0, 3, 0.01)

figuras = lapply(1:6, function(n) 1 - pgompertz(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)


survi = ggplot() +
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C1, colour = "1", linetype = "1"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C2, colour = "2", linetype = "2"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C3, colour = "3", linetype = "3"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C4, colour = "4", linetype = "4"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C5, colour = "5", linetype = "5"),
            size = 1)+
  geom_line(data=dados_risco, 
            mapping=aes(x=x, y=C6, colour = "6", linetype = "6"),
            size = 1)+
  theme_minimal()+
  labs(y = 'S(t)', x = 't') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "tomato",
                       "gold",
                       "royalblue4"),
                     labels = c("a = 3; b = 0,5", 
                                "a = 2; b = 0,5", 
                                "a = 1; b = 0,5",
                                "a = 0,5; b = 0,25",
                                "a = 1; b = 1",
                                "a = 2; b = 2")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("a = 3; b = 0,5", 
                                   "a = 2; b = 0,5", 
                                   "a = 1; b = 0,5",
                                   "a = 0,5; b = 0,25",
                                   "a = 1; b = 1",
                                   "a = 2; b = 2"))

ggsave('figuras/surv_gompertz.pdf', survi, units = 'in', width = 7, height = 5)
