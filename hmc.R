library('rstan')
library('ggplot2')

## Códigos -----

### TGCA ----

dados_tg = read.csv('dados_tg.csv')

#### Modelo usual

fit_usual = stan(file = 'codigos_stan/usu_tgca.stan',
            data = list(N = nrow(dados_tg), T = dados_tg$tempo/365, D = dados_tg$cens), 
            iter = 10000, chains = 1, cores = 1, seed = 154)

print(fit_usual, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

traceplot(fit_usual, pars=c("a", "b"), inc_warmup = TRUE, nrow = 3)

plot(fit_usual)

fit_usual_summary = summary(fit_usual)

fit_usual_summary$summary

fit_usual_summary$summary[,1]

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

ggsave(filename = 'figuras/tgca_parametros_traco.pdf', units = 'in', width = 8, height = 5)

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

ggsave(filename = 'figuras/tgca_parametros_densidade.pdf', units = 'in', width = 8, height = 5)

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

############### Comparando -----

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
    mapping=aes(x=x, y = pgompertz(x, fit_usual_summary$summary[,1]['a'], 
                                   fit_usual_summary$summary[,1]['b'], lower.tail = F),
                colour = "b", linetype = "b"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y= fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'], 
                                                                                                            fit_summary$summary[,1]['b'], lower.tail = F),
                colour = "c", linetype = "c"), 
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, fit_summary_def$summary[,1]['a'], 
                                   fit_summary_def$summary[,1]['b'], lower.tail = F),
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

tgca_tab_resumo = data.frame(a = c(fit_usual_summary$summary[,1]['a'],
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

tgca_tab_resumo |> xtable::xtable(digits = 4) |> print(include.rownames = F)

## Outros conjuntos de dados ----

#### Colon -----
data(cancer, package="survival")

plot(survival::survfit(survival::Surv(time, status) ~ 1, data = colon))

fit = stan( model_code= mix,
            data=list(N = nrow(colon), T = colon$time/365, D = colon$status), 
            iter=5000, chains=2, cores=2)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_def = stan( model_code = defective,
                data=list(N = nrow(colon), T = colon$time/365, D = colon$status), 
                iter=5000, chains=2, cores=2)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)


fit_summary_def = summary(fit_def)

x = colon$time/365

kaplan_meier_s = survival::survfit(survival::Surv(time/365, status) ~ 1, data = colon)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "Kaplan-Meier"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y= fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'], 
                                                                                                                      fit_summary$summary[,1]['b'], lower.tail = F),
                colour = "Com Mix"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, fit_summary_def$summary[,1]['a'], 
                                             fit_summary_def$summary[,1]['b'], lower.tail = F),
                colour = "Defectivo"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom')


#### Retinopatia -----
data(retinopathy, package="survival")

plot(survival::survfit(survival::Surv(futime, status) ~ 1, data = retinopathy))

fit = stan( model_code= mix,
            data=list(N = nrow(retinopathy), T = retinopathy$futime, D = retinopathy$status), 
            iter=5000, chains=2, cores=2)

print(fit, pars=c("a", "b", "theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit)

fit_summary = summary(fit)

fit_summary$summary

fit_def = stan( model_code = defective,
                data=list(N = nrow(retinopathy), T = retinopathy$futime, D = retinopathy$status), 
                iter=5000, chains=2, cores=2)

print(fit_def, pars=c("a", "b"),
      probs=c(0.1, 0.5, 0.9), digits = 3)

plot(fit_def)


fit_summary_def = summary(fit_def)

x = retinopathy$futime

kaplan_meier_s = survival::survfit(survival::Surv(futime, status) ~ 1, data = retinopathy)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() + 
  geom_line(aes(x = Tempo, y = Sobrevivência, colour = "Kaplan-Meier"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
    mapping=aes(x=x, y= fit_summary$summary[,1]['theta'] + (1 - fit_summary$summary[,1]['theta'])*pgompertz(x, fit_summary$summary[,1]['a'], 
                                                                                                                      fit_summary$summary[,1]['b'], lower.tail = F),
                colour = "Com Mix"),
    size = 1) +
  geom_line( 
    mapping=aes(x=x, y = flexsurv::pgompertz(x, fit_summary_def$summary[,1]['a'], 
                                             fit_summary_def$summary[,1]['b'], lower.tail = F),
                colour = "Defectivo"),
    size = 1) +
  ylim(0,1) +
  theme(legend.position = 'bottom')
