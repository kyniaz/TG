library("survival")
library("ggplot2")
source("gompertz.R")

######### Figura ilustrativas ----



############################ S(t) Defectivo

lista_a_bs = list()

lista_a_bs[[1]] = c(-3,2)
lista_a_bs[[2]] = c(-2,1)
lista_a_bs[[3]] = c(-1,1)
lista_a_bs[[4]] = c(-1,2)
lista_a_bs[[5]] = c(1,2)
lista_a_bs[[6]] = c(2,2)


x = seq(0.001, 1, 0.001)

figuras = lapply(1:6, function(n) hazard(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

dados_risco = data.frame(figuras)
colnames(dados_risco) = paste0('C',1:6)

dados_risco = cbind.data.frame(x, dados_risco)

## Sobrevivência ----

x = seq(0, 3, 0.01)

figuras = lapply(1:6, function(n) 1 - flexsurv::pgompertz(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2]))

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
                     labels = c("a = -3; b = 2", 
                                "a = -2; b = 1", 
                                #"1,0.5",
                                "a = -1; b = 1",
                                #"1,0.25",
                                #"1,0.75",
                                "a = -1; b = 2",
                                "a = 1; b = 2",
                                "a = 2; b = 2")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed","twodash",
                                            "longdash", "dotdash"),
                        labels = c("a = -3; b = 2", 
                                   "a = -2; b = 1", 
                                   #"1,0.5",
                                   "a = -1; b = 1",
                                   #"1,0.25",
                                   #"1,0.75",
                                   "a = -1; b = 2",
                                   "a = 1; b = 2",
                                   "a = 2; b = 2"))


ggsave('figuras/surv_gompertz_def.pdf', survi, units = 'in', width = 7, height = 5)


x = seq(0, 5, 0.001)
s = flexsurv::pgompertz(x, shape = -2, rate = 1, lower.tail = F)
H = -log(s)
h = flexsurv::dgompertz(x, shape = -2, rate = 1)/s

dados_plot = data.frame(x, s, H, h)

ggplot() + 
  geom_line(aes(x = x, y = s, colour = "a", linetype = "a"), data = dados_plot, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo', y = " ") +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = h,
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = H,
                colour = "c", linetype = "c"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "tomato"),
                     labels = c("F. de sobrevivência", 
                                "Risco instantâneo", 
                                #"1,0.5",
                                "Risco acumulado")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed"),
                        labels = c("F. de sobrevivência", 
                                   "Risco instantâneo", 
                                   #"1,0.5",
                                   "Risco acumulado"))
  
ggsave(filename = 'figuras/exemplo_def.pdf', units = 'in', width = 7, height = 5)

########################## ESTIMATIVAS FREQUENTISTAS ###################

#### 1 - tgca, 2- colon, 3- melanoma, 4 - retinopathy.

### TGCA ----------
dados_tg = read.csv('dados_tg.csv', header = T)

tempos = dados_tg$tempo/365
cens = dados_tg$cens

npcure::testmz(tempos, cens)

###### Normal ----
log_veros = function(par){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

par_init = c(0.05, 4)

optim_usual = optim(par_init, log_veros, 
                     control = list(fnscale = -1, maxit = 500),
                     method="L-BFGS-B", lower = c(0.0000001, 0.0000001), upper = c(Inf, Inf))

optim_usual

###### Com Mistura ----
log_veros_mix = function(par){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
    (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

par_init = c(0.1, 4, 0.5)

optim_mix = optim(par_init, log_veros_mix, 
                     control = list(fnscale = -1, maxit = 500),
                     method="L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001), upper = c(Inf, Inf, 1))

optim_mix

###### Defeituoso ----
log_veros_def = function(par){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
           (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

par_init = c(-0.05, 1)

optim_def = optim(par_init, log_veros_def, 
                     control = list(fnscale = -1, maxit = 500),
                     method="L-BFGS-B", lower = c(-Inf, 0.0001), upper = c(Inf, Inf))

optim_def

####Df estimativas
estimativas = data.frame(a = c(optim_usual$par[1], optim_mix$par[1], optim_def$par[1]),
                         b = c(optim_usual$par[2], optim_mix$par[2], optim_def$par[2]),
                         p = c('-', optim_mix$par[3], '-'),
                         l = c(optim_usual$value, optim_mix$value, optim_def$value))
                           

xtable::xtable(estimativas, align = 'ccccc', digits = 4)
#veros
# b = 0.031796
# a = seq(-0.2, 0.2, 0.001)
# 
# vals = numeric(length(a))
# 
# for (i in 1:length(a)){
#   vals[i] = log_veros(c(a[i], b))
# }
# 
# plot(a, vals, type = "l")
# 
# plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_tg))

data(cancer, package="survival")

############# Gráficos e Figuras -----


#AIC;
4 + -2*optim_usual$value
4 + -2*optim_def$value
6 + -2*optim_mix$value


#Figuras

x = dados_tg$tempo/365

kaplan_meier_s = survfit(Surv(tempo/365, cens) ~ 1, data = dados_tg)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a",
                ), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x, optim_usual$par[1], optim_usual$par[2], lower.tail = F),
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = optim_mix$par[3] + (1 - optim_mix$par[3])*pgompertz(x, optim_mix$par[1], 
                                                      optim_mix$par[2], lower.tail = F),
                colour = "c", linetype = "c"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),         legend.key.height = unit(0.5, 'cm')) + 
  scale_color_manual(name = "",
                     values = c(
                       "royalblue",
                       "springgreen3",
                       "brown2"),
                     labels = c("Kaplan-Meier","Usual/Defeituoso",
                                "Com Mistura")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted", "dashed"),
                        labels = c("Kaplan-Meier", 
                                   "Usual/Defeituoso", 
                                   "Com Mistura"))


ggsave(filename = 'figuras/tgca_freq.pdf', units = 'in', width = 7, height = 5)

### Colon ----------
data(cancer, package="survival")

tempos = colon$time/365
cens = colon$status

npcure::testmz(tempos, cens)

###### Normal ----
log_veros = function(par){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

par_init = c(1000, 10000)

optim_usual = optim(par_init, log_veros, 
                    control = list(fnscale = -1, maxit = 500),
                    method="L-BFGS-B", lower = c(0.0000001, 0.0000001), upper = c(Inf, Inf))

optim_usual

###### Com Mistura ----
log_veros_mix = function(par){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
               (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

par_init = c(3, 10, 0.2)

optim_mix = optim(par_init, log_veros_mix, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001), upper = c(Inf, Inf, 1))

optim_mix

###### Defeituoso ----
log_veros_def = function(par){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
               (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

par_init = c(-0.05, 1)

optim_def = optim(par_init, log_veros_def, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(-Inf, 0.0001), upper = c(Inf, Inf))

optim_def

#veros

############# Gráficos e Figuras -----


#AIC;
4 + -2*optim_usual$value
4 + -2*optim_def$value
6 + -2*optim_mix$value


#Figuras

x = tempos

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, optim_def$par[1], optim_def$par[2], lower.tail = F),
                colour = "c", linetype = "c"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x, optim_usual$par[1], optim_usual$par[2], lower.tail = F),
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = optim_mix$par[3] + (1 - optim_mix$par[3])*pgompertz(x, optim_mix$par[1], 
                                                                             optim_mix$par[2], lower.tail = F),
                colour = "d", linetype = "d"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),         legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "royalblue",
                       "springgreen3",
                       'orange',
                       "brown2"),
                     labels = c("Kaplan-Meier","Usual", 
                                "Defeituoso",
                                "Com Mistura")) +
  scale_linetype_manual(name = "", values=c("solid", "twodash", 
                                            "dotted", "dashed"),
                        labels = c("Kaplan-Meier", 
                                   "Usual",
                                   "Defeituoso", 
                                   "Com Mistura"))

ggsave(filename = 'figuras/colon_freq.pdf', units = 'in', width = 7, height = 5)

estimativas = data.frame(a = c(optim_usual$par[1], optim_mix$par[1], optim_def$par[1]),
                         b = c(optim_usual$par[2], optim_mix$par[2], optim_def$par[2]),
                         p = c('-', optim_mix$par[3], '-'),
                         l = c(optim_usual$value, optim_mix$value, optim_def$value))


xtable::xtable(estimativas, align = 'ccccc', digits = 4)

### Melanoma -------

melanoma = read.table("melanoma.txt", header = T)

e1690 = melanoma[melanoma$STUDY == '1690' & melanoma$BRESLOW != '.',]

tempos = e1690$SURVTIME
cens = ifelse(e1690$SURVCENS == 2, 1, 0)

npcure::testmz(tempos, cens)

###### Normal ----
log_veros = function(par){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

par_init = c(100, 1000)

optim_usual = optim(par_init, log_veros, 
                    control = list(fnscale = -1, maxit = 500),
                    method="L-BFGS-B", lower = c(0.0000001, 0.0000001), upper = c(Inf, Inf))

optim_usual

###### Com Mistura ----
log_veros_mix = function(par){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
               (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

par_init = c(0.1, 4, 0.5)

optim_mix = optim(par_init, log_veros_mix, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001), upper = c(Inf, Inf, 1))

optim_mix

###### Defeituoso ----
log_veros_def = function(par){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
               (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

par_init = c(-0.05, 1)

optim_def = optim(par_init, log_veros_def, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(-Inf, 0.0001), upper = c(Inf, Inf))

optim_def

############# Gráficos e Figuras -----


#AIC;
4 + -2*optim_usual$value
4 + -2*optim_def$value
6 + -2*optim_mix$value


#Figuras

x = tempos

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, optim_def$par[1], optim_def$par[2], lower.tail = F),
                colour = "c", linetype = "c"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x, optim_usual$par[1], optim_usual$par[2], lower.tail = F),
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = optim_mix$par[3] + (1 - optim_mix$par[3])*pgompertz(x, optim_mix$par[1], 
                                                                             optim_mix$par[2], lower.tail = F),
                colour = "d", linetype = "d"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),         legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "royalblue",
                       "springgreen3",
                       'orange',
                       "brown2"),
                     labels = c("Kaplan-Meier","Usual", 
                                "Defeituoso",
                                "Com Mistura")) +
  scale_linetype_manual(name = "", values=c("solid", "twodash", 
                                            "dotted", "dashed"),
                        labels = c("Kaplan-Meier", 
                                   "Usual",
                                   "Defeituoso", 
                                   "Com Mistura"))

ggsave(filename = 'figuras/melanoma_freq.pdf', units = 'in', width = 7, height = 5)

estimativas = data.frame(a = c(optim_usual$par[1], optim_mix$par[1], optim_def$par[1]),
                         b = c(optim_usual$par[2], optim_mix$par[2], optim_def$par[2]),
                         p = c('-', optim_mix$par[3], '-'),
                         l = c(optim_usual$value, optim_mix$value, optim_def$value))


xtable::xtable(estimativas, align = 'ccccc', digits = 4)

### Diabetic -------

tempos = diabetic$time
cens = diabetic$status

npcure::testmz(tempos, cens)

tempos = flexsurv::bc$rectime/365
cens = flexsurv::bc$censrec
     
###### Normal ----
log_veros = function(par){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

par_init = c(10, 100)

optim_usual = optim(par_init, log_veros, 
                    control = list(fnscale = -1, maxit = 500),
                    method="L-BFGS-B", lower = c(0.0000001, 0.0000001), upper = c(Inf, Inf))

optim_usual

###### Com Mistura ----
log_veros_mix = function(par){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
               (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

par_init = c(0.1, 2, 0.5)

optim_mix = optim(par_init, log_veros_mix, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001), upper = c(Inf, Inf, 1))

optim_mix

###### Defeituoso ----
log_veros_def = function(par){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
               (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

par_init = c(-0.05, 1)

optim_def = optim(par_init, log_veros_def, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(-Inf, 0.0001), upper = c(Inf, Inf))

optim_def

#veros

############# Gráficos e Figuras -----


#AIC;
4 + -2*optim_usual$value
4 + -2*optim_def$value
6 + -2*optim_mix$value


#Figuras

x = tempos

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, optim_def$par[1], optim_def$par[2], lower.tail = F),
                colour = "c", linetype = "c"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x, optim_usual$par[1], optim_usual$par[2], lower.tail = F),
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = optim_mix$par[3] + (1 - optim_mix$par[3])*pgompertz(x, optim_mix$par[1], 
                                                                             optim_mix$par[2], lower.tail = F),
                colour = "d", linetype = "d"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),         legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "royalblue",
                       "springgreen3",
                       'orange',
                       "brown2"),
                     labels = c("Kaplan-Meier","Usual", 
                                "Defeituoso",
                                "Com Mistura")) +
  scale_linetype_manual(name = "", values=c("solid", "twodash", 
                                            "dotted", "dashed"),
                        labels = c("Kaplan-Meier", 
                                   "Usual",
                                   "Defeituoso", 
                                   "Com Mistura"))

ggsave(filename = 'figuras/diabetic_freq.pdf', units = 'in', width = 7, height = 5)

estimativas = data.frame(a = c(optim_usual$par[1], optim_mix$par[1], optim_def$par[1]),
                         b = c(optim_usual$par[2], optim_mix$par[2], optim_def$par[2]),
                         p = c('-', optim_mix$par[3], '-'),
                         l = c(optim_usual$value, optim_mix$value, optim_def$value))


xtable::xtable(estimativas, align = 'ccccc', digits = 4)

### bc -------

bc = flexsurv::bc
# tempos = bc$rectime/365
# cens = bc$censrec

# tempos = stanford2$time/365
# cens = stanford2$status

tempos = ovarian$futim/265
cens = ovarian$fustat

###### Normal ----
log_veros = function(par){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

par_init = c(1, 5)

optim_usual = optim(par_init, log_veros, 
                    control = list(fnscale = -1, maxit = 500),
                    method="L-BFGS-B", lower = c(0.0000001, 0.0000001), upper = c(Inf, Inf))

optim_usual

###### Com Mistura ----
log_veros_mix = function(par){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
               (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

par_init = c(0.05, 2, 0.5)

optim_mix = optim(par_init, log_veros_mix, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(0.0000001, 0.0000001, 0.0000001), upper = c(Inf, Inf, 1))

optim_mix

###### Defeituoso ----
log_veros_def = function(par){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
               (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

par_init = c(-0.05, 1)

optim_def = optim(par_init, log_veros_def, 
                  control = list(fnscale = -1, maxit = 500),
                  method="L-BFGS-B", lower = c(-Inf, 0.0001), upper = c(Inf, Inf))

optim_def

#veros

############# Gráficos e Figuras -----


#AIC;
4 + -2*optim_usual$value
4 + -2*optim_def$value
6 + -2*optim_mix$value


#Figuras

x = tempos

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, optim_def$par[1], optim_def$par[2], lower.tail = F),
                colour = "c", linetype = "c"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x, optim_usual$par[1], optim_usual$par[2], lower.tail = F),
                colour = "b", linetype = "b"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = optim_mix$par[3] + (1 - optim_mix$par[3])*pgompertz(x, optim_mix$par[1], 
                                                                             optim_mix$par[2], lower.tail = F),
                colour = "d", linetype = "d"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom', legend.text = element_text(size=12), legend.key.width= unit(1.5, 'cm'),         legend.key.height = unit(0.5, 'cm')) +
  scale_color_manual(name = "",
                     values = c(
                       "royalblue",
                       "springgreen3",
                       'orange',
                       "brown2"),
                     labels = c("Kaplan-Meier","Usual", 
                                "Defeituoso",
                                "Com Mistura")) +
  scale_linetype_manual(name = "", values=c("solid", "twodash", 
                                            "dotted", "dashed"),
                        labels = c("Kaplan-Meier", 
                                   "Usual",
                                   "Defeituoso", 
                                   "Com Mistura"))

ggsave(filename = 'figuras/ovarian_freq.pdf', units = 'in', width = 7, height = 5)

estimativas = data.frame(a = c(optim_usual$par[1], optim_mix$par[1], optim_def$par[1]),
                         b = c(optim_usual$par[2], optim_mix$par[2], optim_def$par[2]),
                         p = c('-', optim_mix$par[3], '-'),
                         l = c(optim_usual$value, optim_mix$value, optim_def$value))


xtable::xtable(estimativas, align = 'ccccc', digits = 4)
