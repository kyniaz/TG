library('rstan')
library('ggplot2')

## Códigos -----
log_veros = function(par, tempos, cens){
  return(sum(cens*dgompertz(tempos, par[1], par[2], ln = T) +
               (1 - cens)*pgompertz(tempos, par[1], par[2], ln = T, lower.tail = F)))
}

log_veros_mix = function(par, tempos, cens){
  return(sum(cens*log(1-par[3]) + cens*dgompertz(tempos, par[1], par[2], ln = T) + 
               (1- cens)*log(par[3]  + (1-par[3])*(1 - pgompertz(tempos, a = par[1], b = par[2])))))
}

log_veros_def = function(par, tempos, cens){
  return(sum(cens*flexsurv::dgompertz(tempos, par[1], par[2], log = T) +
               (1 - cens)*flexsurv::pgompertz(tempos, par[1], par[2], log = T, lower.tail = F)))
}

calcula_dic = function(tempos, cens, log_veros, cadeia_a, cadeia_b, cadeia_p = NULL){
  a = mean(cadeia_a)
  b = mean(cadeia_b)
  
  if(!is.null(cadeia_p)){
    p = mean(cadeia_p)
    par = c(a,b,p)
    aux = numeric(length = length(cadeia_a))
    
    for(i in 1:length(cadeia_a)){
      aux[i] = log_veros(c(cadeia_a[i],cadeia_b[i], cadeia_p[i]), tempos, cens)
    }
    
    p_dic = 2*(log_veros(par, tempos, cens) - mean(aux))
    
    dic = -2*log_veros(par, tempos, cens) + 2*p_dic
  }
  else{
    par = c(a,b)
    aux = numeric(length = length(cadeia_a))
    
    for(i in 1:length(cadeia_a)){
      aux[i] = log_veros(c(cadeia_a[i],cadeia_b[i]), tempos, cens)
    }
    
    p_dic = 2*(log_veros(par, tempos, cens) - mean(aux))
    
    dic = -2*log_veros(par, tempos, cens) + 2*p_dic
  }
  return(dic)
}
            
### TGCA ----

dados_tg = read.csv('dados_tg.csv')

#### Modelo usual

fit_usual = stan(file = 'codigos_stan/usu_tgca.stan',
            data = list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
            iter = 11000, warmup = 1000, chains = 1, cores = 1, seed = 154)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), window = c(1,2000), inc_warmup  = T, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

# if(length(cadeias) > 1){
#   wup = attr(cadeias[[1]],"args")$warmup
#   ite = attr(cadeias[[1]],"args")$iter
#   cadeia_conjunta = numeric((wup - ite)*length(cadeias))
#   
#   for(i in 1:length(cadeias)){
#     indices = seq(i, (wup - ite)*length(cadeias), length(cadeias))
#     cadeia_conjunta[indices] = cadeias[[1]]
#   }
# } else {
#   cadeias = fit_usual@sim$samples
# }

tgca_a = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
tgca_b = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(dados_tg$tempo/365, dados_tg$cens, log_veros,
           tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
           tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |>
  print(digits = 22)

#1222.423

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b)
#### Modelo com Mistura ----

fit = stan( file = 'codigos_stan/mix_tgca.stan',
            data = list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit, pars=c("a", "b", "theta"),                                                   
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit, pars=c("a", "b", "theta"), inc_warmup = TRUE, nrow = 3)

plot(fit)

cadeias = fit@sim$samples

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b, p = cadeias[[1]]$theta)

tgca_a2 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
tgca_b2 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))
tgca_p2 = median(tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

calcula_dic(dados_tg$tempo/365, dados_tg$cens, log_veros_mix, 
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2)) |>
  print(digits = 22)


fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

##### Modelo Defectivo -----

fit_def = stan(file = 'codigos_stan/def_tgca.stan',
            data=list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
            iter=10000, chains=1, cores=1)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)

fit_summary_def = summary(fit_def)

fit_summary_def$summary

fit_summary_def$summary[,1]

cadeias = fit_def@sim$samples

tgca_a3 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
tgca_b3 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(dados_tg$tempo/365, dados_tg$cens, log_veros_def,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |>
  print(digits = 22)

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b)

library("survival")

x = dados_tg$tempo/365

kaplan_meier_s = survfit(Surv(tempo/365, cens) ~ 1, data = dados_tg)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), 
            data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x,tgca_a , #fit_usual_summary$summary[,1]['a']
                                   tgca_b, lower.tail = F), #fit_usual_summary$summary[,1]['b']
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, tgca_a3, #fit_summary_def$summary[,1]['a']
                                   tgca_b3, lower.tail = F), #fit_summary_def$summary[,1]['b']
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= tgca_p2 + (1 - tgca_p2)*pgompertz(x, tgca_a2, tgca_b2, lower.tail = F), #fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'],     #fit_summary$summary[,1]['b']
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

ggsave(filename = 'figuras/tgca_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC

tgca_tab_resumo = data.frame(a = c(tgca_a, #fit_usual_summary$summary[,1]['a'],
                                   tgca_a2,#fit_summary$summary[,1]['a'],
                                   tgca_a3), #fit_summary_def$summary[,1]['a']),
                             b = c(tgca_b, #fit_usual_summary$summary[,1]['b'],
                                   tgca_b2,#fit_summary$summary[,1]['b'],
                                   tgca_b3),#fit_summary_def$summary[,1]['b']),
                             p = c('',
                               tgca_p2,#fit_summary$summary[,1]['theta'],
                               ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))

tgca_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)


## Outros conjuntos de dados ----

#### Colon -----
data(cancer, package="survival")

fit_usual = stan(file = 'codigos_stan/usu_colon.stan',
                 data = list(N = nrow(colon), T = colon$time/365, D = colon$status), 
                 iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), inc_warmup = TRUE, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

cadeias = fit_usual@sim$samples

colon_a = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
colon_b = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(colon$time/365, colon$status, log_veros,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL)

#5786.717

#### Modelo com Mistura ----

fit = stan( file = 'codigos_stan/mix_colon.stan',
            data = list(N = nrow(colon), T = colon$time/365, D = colon$status),  
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit, pars=c("a", "b", "theta"), inc_warmup = TRUE, nrow = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

cadeias = fit@sim$samples

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b, p = cadeias[[1]]$theta)

colon_a2 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
colon_b2 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))
colon_p2 = median(tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

calcula_dic(colon$time/365, colon$status, log_veros_mix, 
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

#5575.265

##### Modelo Defectivo -----

fit_def = stan(file = 'codigos_stan/def_colon.stan',
               data = list(N = nrow(colon), T = colon$time/365, D = colon$status),  
               iter=10000, chains=1, cores=1)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)

fit_summary_def = summary(fit_def)

fit_summary_def$summary

fit_summary_def$summary[,1]

cadeias = fit_def@sim$samples

colon_a3 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
colon_b3 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(colon$time/365, colon$status, log_veros_def,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |>
  print(digits = 22)

# 5585.222
############### Comparando -----

library("survival")

x = colon$time/365

kaplan_meier_s = survfit(Surv(time/365, status) ~ 1, data = colon)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), 
            data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x,colon_a , #fit_usual_summary$summary[,1]['a']
                                   colon_b, lower.tail = F), #fit_usual_summary$summary[,1]['b']
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, colon_a3, #fit_summary_def$summary[,1]['a']
                                             colon_b3, lower.tail = F), #fit_summary_def$summary[,1]['b']
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= colon_p2 + (1 - colon_p2)*pgompertz(x, colon_a2, colon_b2, lower.tail = F), #fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'],     #fit_summary$summary[,1]['b']
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

ggsave(filename = 'figuras/colon_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC

colon_tab_resumo = data.frame(a = c(colon_a, #fit_usual_summary$summary[,1]['a'],
                                   colon_a2,#fit_summary$summary[,1]['a'],
                                   colon_a3), #fit_summary_def$summary[,1]['a']),
                             b = c(colon_b, #fit_usual_summary$summary[,1]['b'],
                                   colon_b2,#fit_summary$summary[,1]['b'],
                                   colon_b3),#fit_summary_def$summary[,1]['b']),
                             p = c('',
                                   colon_p2,#fit_summary$summary[,1]['theta'],
                                   ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))

colon_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)


### Diabetic  -----
data(cancer, package="survival")

fit_usual = stan(file = 'codigos_stan/usu_diab.stan',
                 data = list(N = nrow(diabetic), T = diabetic$time/12, D = diabetic$status), 
                 iter = 10000, chains = 1, cores = 1, seed = 245)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), inc_warmup = TRUE, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

cadeias = fit_usual@sim$samples

dia_a = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
dia_b = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(diabetic$time/12, diabetic$status, log_veros,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |>
  print(digits = 22)

#936.3835

##### Modelo com Mistura ----

fit = stan( file = 'codigos_stan/mix_diab.stan',
            data = list(N = nrow(diabetic), T = diabetic$time/12, D = diabetic$status), 
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit, pars=c("a", "b", "theta"), inc_warmup = TRUE, nrow = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

cadeias = fit@sim$samples

dia_a2 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
dia_b2 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))
dia_p2 = median(tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b, p = cadeias[[1]]$theta)

calcula_dic(diabetic$time/12, diabetic$status, log_veros_mix, 
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2)) |> print(digits = 22)
#922.2386

##### Modelo Defectivo -----

fit_def = stan(file = 'codigos_stan/def_diab.stan',
               data = list(N = nrow(diabetic), T = diabetic$time/12, D = diabetic$status), 
               iter=10000, chains=1, cores=1)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)

fit_summary_def = summary(fit_def)

fit_summary_def$summary

fit_summary_def$summary[,1]

cadeias = fit_def@sim$samples

dia_a3 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
dia_b3 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(diabetic$time/12, diabetic$status, log_veros_def,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |> print(digits = 22)

#  922.3669
############### Comparando -----

library("survival")

x = diabetic$time/12

kaplan_meier_s = survfit(Surv(time/12, status) ~ 1, data = diabetic)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), 
            data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x,dia_a , #fit_usual_summary$summary[,1]['a']
                                   dia_b, lower.tail = F), #fit_usual_summary$summary[,1]['b']
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, dia_a3, #fit_summary_def$summary[,1]['a']
                                             dia_b3, lower.tail = F), #fit_summary_def$summary[,1]['b']
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= dia_p2 + (1 - dia_p2)*pgompertz(x, dia_a2, dia_b2, lower.tail = F), #fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'],     #fit_summary$summary[,1]['b']
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

ggsave(filename = 'figuras/dia_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC

dia_tab_resumo = data.frame(a = c(dia_a, #fit_usual_summary$summary[,1]['a'],
                                   dia_a2,#fit_summary$summary[,1]['a'],
                                   dia_a3), #fit_summary_def$summary[,1]['a']),
                             b = c(dia_b, #fit_usual_summary$summary[,1]['b'],
                                   dia_b2,#fit_summary$summary[,1]['b'],
                                   dia_b3),#fit_summary_def$summary[,1]['b']),
                             p = c('',
                                   dia_p2,#fit_summary$summary[,1]['theta'],
                                   ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))

ggsave(filename = 'figuras/diabetic_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC


diab_tab_resumo = data.frame(a = c(fit_usual_summary$summary[,1]['a'],
                                   fit_summary$summary[,1]['a'],
                                   fit_summary_def$summary[,1]['a']),
                             b = c(fit_usual_summary$summary[,1]['b'],
                                   fit_summary$summary[,1]['b'],
                                   fit_summary_def$summary[,1]['b']),
                             p = c('',
                                   fit_summary$summary[,1]['theta'],
                                   ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))

diab_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)

### Diabetic


### Ovarian ----

#### Modelo usual

fit_usual = stan(file = 'codigos_stan/usu_ova.stan',
                 data = list(N = nrow(ovarian), T = ovarian$futime/365, D = ovarian$fustat), 
                 iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), inc_warmup = TRUE, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

cadeias = fit_usual@sim$samples

ova_a = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
ova_b = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(ovarian$futime/365, ovarian$fustat, log_veros,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL)

#55.69841

#### Modelo com Mistura ----

fit = stan( file = 'codigos_stan/mix_ova.stan',
            data = list(N = nrow(ovarian), T = ovarian$futime/365, D = ovarian$fustat), 
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit, pars=c("a", "b", "theta"), inc_warmup = TRUE, nrow = 3)

plot(fit)

cadeias = fit@sim$samples

ova_a2 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
ova_b2 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))
ova_p2 = median(tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b, p = cadeias[[1]]$theta)

calcula_dic(ovarian$futime/365, ovarian$fustat, log_veros_mix, 
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

#52.68733
###### Analise de convergência ----

a = ggplot() + 
  geom_line(aes(x = index, y = a), data = cadeias_df,
            col = 'royalblue') +
  labs(x = '', y = expression(hat(a))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

b = ggplot() + 
  geom_line(aes(x = index, y = b), data = cadeias_df,
            col = 'springgreen3') +
  labs(x = '', y = expression(hat(b))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

p = ggplot() + 
  geom_line(aes(x = index, y = p), data = cadeias_df,
            col = 'tomato') +
  labs(x = 'Ordem', y = expression(hat(pi))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

abp_traco = gridExtra::grid.arrange(a,b,p, nrow = 3)

ggsave(filename = 'figuras/ova_parametros_traco.pdf', units = 'in', width = 8, height = 5)

###### Densidades ------

a = ggplot() + 
  geom_histogram(aes(x = a, y = ..density..), data = cadeias_df, bins = 30,
                 fill = 'royalblue', col = 'grey') +
  labs(x = expression(hat(a)), y = 'Densidade') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

b = ggplot() + 
  geom_histogram(aes(x = b, y = ..density..), data = cadeias_df, bins = 30,
                 fill = 'springgreen3', col = 'grey') +
  labs(x = expression(hat(b)), y = '') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

p = ggplot() + 
  geom_histogram(aes(x = p, y = ..density..), data = cadeias_df, bins = 30,
                 fill = 'tomato', col = 'grey') +
  labs(x = expression(hat(pi)), y = '') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) 

abp = gridExtra::grid.arrange(a,b,p, ncol = 3)

ggsave(filename = 'figuras/ova_parametros_densidade.pdf', units = 'in', width = 8, height = 5)

fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

##### Modelo Defectivo -----

fit_def = stan(file = 'codigos_stan/def_ova.stan',
               data = list(N = nrow(ovarian), T = ovarian$futime/365, D = ovarian$fustat), 
               iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)

fit_summary_def = summary(fit_def)

fit_summary_def$summary

fit_summary_def$summary[,1]

cadeias = fit_def@sim$samples

ova_a3 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
ova_b3 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(ovarian$futime/365, ovarian$fustat, log_veros_def,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL)

#57.80024

############### Comparando -----

library("survival")

x = ovarian$futime/365

kaplan_meier_s = survfit(Surv(futime/365, fustat) ~ 1, data = ovarian)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), 
            data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x,ova_a , #fit_usual_summary$summary[,1]['a']
                                   ova_b, lower.tail = F), #fit_usual_summary$summary[,1]['b']
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, ova_a3, #fit_summary_def$summary[,1]['a']
                                             ova_b3, lower.tail = F), #fit_summary_def$summary[,1]['b']
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= ova_p2 + (1 - ova_p2)*pgompertz(x, ova_a2, ova_b2, lower.tail = F), #fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'],     #fit_summary$summary[,1]['b']
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

ggsave(filename = 'figuras/ova_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC

ova_tab_resumo = data.frame(a = c(ova_a, #fit_usual_summary$summary[,1]['a'],
                                   ova_a2,#fit_summary$summary[,1]['a'],
                                   ova_a3), #fit_summary_def$summary[,1]['a']),
                             b = c(ova_b, #fit_usual_summary$summary[,1]['b'],
                                   ova_b2,#fit_summary$summary[,1]['b'],
                                   ova_b3),#fit_summary_def$summary[,1]['b']),
                             p = c('',
                                   ova_p2,#fit_summary$summary[,1]['theta'],
                                   ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))


ova_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)


### Melanoma  -----
melanoma = read.table("melanoma.txt", header = T)

e1690 = melanoma[melanoma$STUDY == '1690' & melanoma$BRESLOW != '.',]

tempos = e1690$SURVTIME
cens = ifelse(e1690$SURVCENS == 2, 1, 0)

fit_usual = stan(file = 'codigos_stan/usu_mel.stan',
                 data = list(N = nrow(e1690), T = tempos, D = cens), 
                 iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), inc_warmup = TRUE, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

cadeias = fit_usual@sim$samples

mel_a = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
mel_b = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(tempos, cens, log_veros,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |>
  print(digits = 22)

#1100.2292

##### Modelo com Mistura ----

fit = stan( file = 'codigos_stan/mix_mel.stan',
            data = list(N = nrow(e1690), T = tempos, D = cens), 
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit, pars=c("a", "b", "theta"), inc_warmup = TRUE, nrow = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_summary$summary[,1]

cadeias = fit@sim$samples

mel_a2 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
mel_b2 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))
mel_p2 = median(tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2))

cadeias_df = data.frame(index = 1:length(cadeias[[1]]$a), a = cadeias[[1]]$a, b = cadeias[[1]]$b, p = cadeias[[1]]$theta)

calcula_dic(tempos, cens, log_veros_mix, 
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), 
            tail(cadeias[[1]]$theta, n = length(cadeias[[1]]$a)/2)) |> print(digits = 22)
#1084.8475

##### Modelo Defectivo -----

fit_def = stan(file = 'codigos_stan/def_mel.stan',
               data = list(N = nrow(e1690), T = tempos, D = cens), 
               iter=10000, chains=1, cores=1)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)

fit_summary_def = summary(fit_def)

fit_summary_def$summary

fit_summary_def$summary[,1]

cadeias = fit_def@sim$samples

mel_a3 = median(tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2))
mel_b3 = median(tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2))

calcula_dic(tempos, cens, log_veros_def,
            tail(cadeias[[1]]$a, n = length(cadeias[[1]]$a)/2),
            tail(cadeias[[1]]$b, n = length(cadeias[[1]]$a)/2), NULL) |> print(digits = 22)

# 1096.4392
############### Comparando -----

library("survival")

x = tempos

kaplan_meier_s = survfit(Surv(tempos, cens) ~ 1)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')


ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a", linetype = "a"), 
            data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y = pgompertz(x,mel_a , #fit_usual_summary$summary[,1]['a']
                                   mel_b, lower.tail = F), #fit_usual_summary$summary[,1]['b']
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, mel_a3, #fit_summary_def$summary[,1]['a']
                                             mel_b3, lower.tail = F), #fit_summary_def$summary[,1]['b']
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= mel_p2 + (1 - mel_p2)*pgompertz(x, mel_a2, mel_b2, lower.tail = F), #fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'],     #fit_summary$summary[,1]['b']
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

ggsave(filename = 'figuras/mel_bayes.pdf', units = 'in', width = 7, height = 5)

#### Comparativo via tabela.

###DIC

mel_tab_resumo = data.frame(a = c(mel_a, #fit_usual_summary$summary[,1]['a'],
                                   mel_a2,#fit_summary$summary[,1]['a'],
                                   mel_a3), #fit_summary_def$summary[,1]['a']),
                             b = c(mel_b, #fit_usual_summary$summary[,1]['b'],
                                   mel_b2,#fit_summary$summary[,1]['b'],
                                   mel_b3),#fit_summary_def$summary[,1]['b']),
                             p = c('',
                                   mel_p2,#fit_summary$summary[,1]['theta'],
                                   ''),
                             l = c(fit_usual_summary$summary[,1]['lp__'],
                                   fit_summary$summary[,1]['lp__'],
                                   fit_summary_def$summary[,1]['lp__']))

mel_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)
