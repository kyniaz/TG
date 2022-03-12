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
  
  start = c(0.1,0.1)
  
  param = optim(start, log_like, control = list(fnscale = -1, maxit = 500),
                method="L-BFGS-B", lower = c(0.001,0.001), upper = c(Inf,Inf))
  
  return(param$par)
}


####### Sem censura ----


set.seed(154)

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5),
  a_b = numeric(5),
  b_b = numeric(5)
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
  
  B = 200
  BB = 500
  a = rep(NA, BB)
  b = rep(NA, BB)
  
  a[1] = 0.5
  b[1] = 0.5
  
  stime = Sys.time()
  for (ii in 2:BB){
    delta = trunc(runif(1, 1, 3))
    a[ii] = ifelse(delta == 1, metropolis_gibbs(a[ii-1], b[ii-1], rprop_gamma, ldprop_gamma, tempos, d , B, start = 0.5), a[ii-1])
    b[ii] = ifelse(delta == 2, metropolis_gibbs(a[ii-1], b[ii-1], rprop_gamma, ldprop_gamma, tempos, d , B, start = 0.5), b[ii-1])
  }
  etime = Sys.time() - stime
  print(paste("Exec time: ", etime))
  
  estimativas[i,c(4,5)] = c(median(a), median(b))
}

estimativas


#Bayesiano

metropolis_gibbs <- function(a, b, rprop, ldprop, tempos, d, B = 10^3, start = 0.5)
{
  chain <- as.list(rep(NA, B))
  chain[[1]] <- start
  
  if(delta == 1) {
    for(ii in 2:B)
    {
      prop <- rprop(chain[[ii-1]])
      
      lratio <- ldtgt1(tempos, d, prop, b) -ldtgt1(tempos, d, chain[[ii-1]], b) +
        ldprop(prop,chain[[ii-1]])-
        ldprop(chain[[ii-1]],prop)
      
      if(log(runif(1)) <= lratio) chain[[ii]] <- prop
      else{
        chain[[ii]] <- chain[[ii-1]]
      }
    }
  }
  if(delta == 2) {
    for(ii in 2:B)
    {
      prop <- rprop(chain[[ii-1]])
      
      lratio <- ldtgt2(tempos, d, a, prop) - ldtgt2(tempos, d, a, chain[[ii-1]]) +
        ldprop(prop,chain[[ii-1]])-
        ldprop(chain[[ii-1]],prop)
      
      if(is.nan(lratio)) chain[[ii]] <- chain[[ii-1]]
      else if(log(runif(1)) <= lratio) chain[[ii]] <- prop
      else{
        chain[[ii]] <- chain[[ii-1]]
      }
    }
  }
  return(tail(unlist(chain), n = 1))
}

######## Logs das Densidades ################

#### com mistura;

gamma_hiper = c(1,0.5)

log_veros = function(tempos, cens, alfa, beta_p){
  return(sum(cens*dgompertz(tempos, a = alfa, b = beta_p, log_opt = T) + 
               (1-cens)*(sgompertz(tempos, a = alfa, b = beta_p,log_opt = T))))
}


ldtgt1 <- function(tempos, d, a, b) {
  #if(prob>1) {print(paste0('delta(1) :', delta, ' e prob: ', prob))}
  return(log_veros(tempos, d, a, b) + log(dgamma(a,gamma_hiper[1],gamma_hiper[2])) + log(dgamma(b,gamma_hiper[1],gamma_hiper[2])))
}

ldtgt2 <- function(tempos, d, a, b) {
  return(log_veros(tempos, d, a, b) + log(dgamma(a,gamma_hiper[1],gamma_hiper[2])) + log(dgamma(b,gamma_hiper[1],gamma_hiper[2])))
}

####### Propostas ----

rprop_gamma = function(ant) rexp(1, 1/ant)
ldprop_gamma = function(ant, prop) dexp(prop, 1/ant, log = T)

grafs <-lapply(1:length(lista_tempos),
               function(col) {
                 
                 dados_simulados = data.frame(tempos = lista_tempos[[col]], cens = rep(1, length(lista_tempos[[col]])))
                 kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)
                 
                 plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
                 colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')
               
                 
                 ggplot2::ggplot(data.frame(x = lista_tempos[[col]]), 
                 aes(x)) +
                 geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'brown2') +
                 #geom_histogram(aes(y=..density..), color="grey", fill="springgreen3", bins = 10) +
                 theme_minimal() +
                 theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom') + 
                 labs(x = 'Tempos', y = 'Sobrevivência') +
                 geom_line(aes(x=x, y = y, colour = "EMV"), data = data.frame(x = lista_tempos[[col]],
                                                                            y = sgompertz(lista_tempos[[col]], 
                                                                            a = estimativas$a[col],
                                                                            b = estimativas$b[col]))) +
                 geom_line(aes(x=x, y = y, colour = "Bayesiana"), data = data.frame(x = lista_tempos[[col]],
                                                                                    y = sgompertz(lista_tempos[[col]], 
                                                                                    a = estimativas$a_b[col],
                                                                                    b = estimativas$b_b[col]))) +
                 scale_color_manual(name = "Método de Estimação",
                                    values = c(
                                      "dodgerblue",
                                      "springgreen3"),
                                    #labels = c("EMV", "Bayesiana")
                                    ) 
               }
              )

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legenda = g_legend(grafs[[1]])

gridExtra::grid.arrange(
  arrangeGrob(grafs[[1]] + theme(legend.position="none"), 
              grafs[[2]] + theme(legend.position="none"), 
              grafs[[3]] + theme(legend.position="none"), 
              grafs[[4]] + theme(legend.position="none"), 
              grafs[[5]] + theme(legend.position="none"), nrow=1),
              legenda, nrow=2, heights=c(10, 1))
#obj

xtable::xtable(estimativas, digits = 3)

dados_simulados = data.frame(tempos, rep(1, max(n)))
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

plot_sim = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'brown2') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #geom_point(aes(x = Tempo, y = Sobrevivência), data = plot_dados[plot_dados$Evento == 1,], shape = 3, color = 'royalblue', size = 1.25) +
  labs(x = 'Tempo') +
  ylim(0, 1) + 
  theme_minimal() 

plot_sim

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

freq_g = plot_sim + geom_line(aes(x=x, y=d, colour = "a"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[2]]) +
  geom_line(aes(x=x, y=d, colour = "c"), data = camadas[[3]]) +
  geom_line(aes(x=x, y=d, colour = "d"), data = camadas[[4]]) +
  geom_line(aes(x=x, y=d, colour = "e"), data = camadas[[5]]) +
  scale_color_manual(name = "Tamanho Amostral",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "darkseagreen2",
                       "gold"),
                     labels = c("25","50","100","200","500")) +
  theme(legend.position = 'bottom')

gridExtra::grid.arrange(freq_g, baye_g, ncol = 2)

hist(tempos, breaks = 20, probability = T)
lines(x_teste, dgompertz(x_teste, mean(a), mean(b)))

# Com censura ----

x_teste = seq(0.001, 6, 0.001)

plot(x_teste, dgompertz(x_teste, 1, 1.35), type = "l", col = "blue")
lines(x_teste, dgompertz(x_teste, 1, 2), col = "red")

n = c(25, 50, 100, 200, 500)

estimativas = data.frame(
  n = numeric(5),
  a = numeric(5),
  b = numeric(5),
  a_b = numeric(5),
  b_b = numeric(5),
  cens = numeric(5)
)

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
  
  B = 200
  BB = 500
  a = rep(NA, BB)
  b = rep(NA, BB)
  
  a[1] = 0.5
  b[1] = 1
  
  stime = Sys.time()
  for (ii in 2:BB){
    delta = trunc(runif(1, 1, 3))
    a[ii] = ifelse(delta == 1, metropolis_gibbs(a[ii-1], b[ii-1], rprop_gamma, 
                                                ldprop_gamma, tempos, eventos , B, start = 0.5), a[ii-1])
    b[ii] = ifelse(delta == 2, metropolis_gibbs(a[ii-1], b[ii-1], rprop_gamma, 
                                                ldprop_gamma, tempos, eventos , B, start = 0.5), b[ii-1])
  }
  
  etime = Sys.time() - stime
  print(paste("Exec time: ", etime))
  
  estimativas[i,c(4,5)] = c(median(a), median(b))
  
  estimativas[i, 6] = 1 - sum(eventos)/n[i]
}

estimativas

dados_simulados = data.frame(tempos, eventos)
colnames(dados_simulados) = c("tempos", "cens")

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1, data = dados_simulados)

plot_dados = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

plot_sim = ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'brown2') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #geom_point(aes(x = Tempo, y = Sobrevivência), data = plot_dados[plot_dados$Evento == 1,], shape = 3, color = 'royalblue', size = 1.25) +
  labs(x = 'Tempo') +
  ylim(0, 1) + 
  theme_minimal() 

plot_sim

estimativas


######Para cada conjunto, plotar a curva encima do km;

xplot = seq(0.01, 3, 0.01)
camadas = list()

for(i in 1:5){
  nova_camada = data.frame(xplot, 1 - pgompertz(xplot, estimativas$a_b[i], 
                                                estimativas$b_b[i]))
  colnames(nova_camada) = c('x','d')
  camadas[[i]] = nova_camada
}

baye_cc = plot_sim + geom_line(aes(x=x, y=d, colour = "a"), data = camadas[[1]]) + 
  geom_line(aes(x=x, y=d, colour = "b"), data = camadas[[2]]) +
  geom_line(aes(x=x, y=d, colour = "c"), data = camadas[[3]]) +
  geom_line(aes(x=x, y=d, colour = "d"), data = camadas[[4]]) +
  geom_line(aes(x=x, y=d, colour = "e"), data = camadas[[5]]) +
  scale_color_manual(name = "Tamanho Amostral",
                     values = c(
                       "dodgerblue",
                       "springgreen3",
                       "lightsteelblue",
                       "darkseagreen2",
                       "gold"),
                     labels = c("25","50","100","200","500")) +
  theme(legend.position = 'bottom')

gridExtra::grid.arrange(baye_cc, freq_cc, ncol = 2)

print(xtable(estimativas, digits = 3),include.rownames = F)

# Risco, Dens e Sobr ----

hazard = function(x, a, b){
  return(dgompertz(x,a,b, log_opt = T) - (sgompertz(x,a,b, log_opt = T)))
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

x = seq(0.001, 3, 0.001)

figuras = lapply(1:6, function(n) exp(hazard(x, lista_a_bs[[n]][1], lista_a_bs[[n]][2])))

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

ggsave('figuras_TG/haz_gompertz.pdf', haz, units = 'in', width = 7, height = 5)
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
