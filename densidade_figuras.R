# Risco, Dens e Sobr ----
source('gompertz.R')
library('ggplot2')

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
