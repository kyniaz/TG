library('survival')

#survcure = function()
dados_tg = read.csv('dados_tg.csv', header = T)

form = as.formula(Surv(tempo/365, cens) ~ 1)
dados = dados_tg

times = dados[[match(all.vars(form)[1], colnames(dados))]]/365
cens = dados[[match(all.vars(form)[2], colnames(dados))]]

#build the X matrix
X = model.matrix(form, dados)

#Optimize
#build the log-loglikelihood

qtde_par = ncol(X)
par_strings_pi = paste0('par[',2:(qtde_par+1),']','*X[[',1:qtde_par,']]')
par_strings_b = paste0('par[',(qtde_par+2):(2*qtde_par+1),']','*X[[',1:qtde_par,']]')

X = data.frame(X)

#Com mistura
log_likelihood = eval(parse(text = paste0('function(par) {
                b_par =  exp(',paste0(par_strings_b, collapse = ' + ',sep = ''),')
                eta = ',paste0(par_strings_pi, collapse = ' + ',sep = ''),'
                mu = 1/(1 + exp(-eta))
                #mu = pnorm(eta)
                l = cens*log(1-mu) + cens*dgompertz(times, par[1], b_par, T) + 
                (1- cens)*log(mu  + (1-mu)*(1 - pgompertz(times, a = par[1], b = b_par)))
                return(sum(l, na.rm = T))
                }')))

#optimize
par_init = numeric(2*ncol(X) + 1)

#par_init = numeric(ncol(X) + 2)

#Set some initial pars
par_init[1] = 0.1
par_init[qtde_par + 2] = 0.5
par_init[2] = 0.5

lower_optim = c(0.00001, rep(-Inf, 2*ncol(X)))
upper_optim = rep(Inf, 2*ncol(X) + 1)

#lower_optim = c(0.00001, rep(-Inf, ncol(X)))
#upper_optim = rep(Inf, ncol(X) + 1)

optim_result = optim(par_init, log_likelihood, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = lower_optim,
              upper = upper_optim)
mix = optim_result$par
aic_mix = -2*optim_result$value + 2*3

#Sem mistura
log_likelihood = eval(parse(text = paste0('function(par) {
                b_par =  exp(',paste0(par_strings_b, collapse = ' + ',sep = ''),')
                eta = ',paste0(par_strings_pi, collapse = ' + ',sep = ''),'
                mu = 1/(1 + exp(-eta))
                #mu = pnorm(eta)
                l = cens*log(-log(mu)) + cens*dgompertz(times, par[1], b_par, ln = T) + 
                (log(mu)*(pgompertz(times, a = par[1], b = b_par)))
                
                return(sum(l, na.rm = T))
                }')))

#optimize
par_init = numeric(2*ncol(X) + 1)

#par_init = numeric(ncol(X) + 2)

#Set some initial pars
par_init[1] = 0.1
par_init[qtde_par + 2] = 0.5
par_init[2] = 0.5

lower_optim = c(0.00001, rep(-Inf, 2*ncol(X)))
upper_optim = rep(Inf, 2*ncol(X) + 1)

#lower_optim = c(0.00001, rep(-Inf, ncol(X)))
#upper_optim = rep(Inf, ncol(X) + 1)

optim_result = optim(par_init, log_likelihood, control = list(fnscale = -1, maxit = 500),
                     method="L-BFGS-B", lower = lower_optim,
                     upper = upper_optim)
optim_result
no_mix = optim_result$par
aic_no_mix = -2*optim_result$value + 2*3

library('jskm')
library('ggplot2')


######### Figuras sem covariáveis
logit_mu = function(eta) return(1/(1+exp(-eta)))

x = dados$tempo/365

kaplan_meier_s = survfit(Surv(tempo/365, cens) ~ 1, data = dados_tg)

dados_km = data.frame(kaplan_meier_s$time, kaplan_meier_s$surv, kaplan_meier_s$n.event)
colnames(dados_km) = c('Tempo', 'Sobrevivência', 'Evento')

#Sem cura
sem_cura = survfit_gompertz(dados$tempo/365, dados$cens)
sem_cura

ggplot() + 
geom_line(aes(x = Tempo, y = Sobrevivência, colour = "a"), data = dados_km, size = 1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Tempo') +
  theme_minimal() +
  geom_line( 
          mapping=aes(x=x, y= logit_mu(mix[2]) + 
                     (1 - logit_mu(mix[2]))*pgompertz(x,mix[1], exp(mix[3]), lower.tail = F), colour = "b"),
          size = 0.5) +
  geom_line( 
            mapping=aes(x=x, y= logit_mu(no_mix[2])^pgompertz(x, no_mix[1], exp(no_mix[3])),
                        colour = "c"),
            size = 0.5) +
  geom_line( 
            mapping=aes(x=x, y= pgompertz(x, sem_cura[1], sem_cura[2], lower.tail = F),
                        colour = "d"),
            size = 0.5) +
  ylim(0,1) +
  scale_color_manual(name = "",
                     values = c(
                       "black",
                       "springgreen3",
                       "royalblue",
                       "brown2"),
                     labels = c("Kaplan-Meier","Com Mistura",
                                "Sem Mistura", "Sem Cura")) +
  theme(legend.position = 'bottom')

ggsave('figuras/curvas_tgca.pdf', units = 'in', height = 5, width = 7)

##AIC


#IDADE
dados_plot = data.frame(x = dados$tempo/365,
                        C1 = logit_mu(optim_result$par[2]) + 
                             (1 - logit_mu(optim_result$par[2]))*
                             sgompertz(dados$tempo/365, optim_result$par[1],
                             exp(optim_result$par[5])),
                        C2 = logit_mu(optim_result$par[2] + optim_result$par[3]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[3]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[5] + optim_result$par[6])),
                        C3 = logit_mu(optim_result$par[2] + optim_result$par[4]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[4]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[5] + optim_result$par[7]))
                        )

jskm(survfit(Surv(tempo/365, cens) ~ trat, data = dados_tg)) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C1, colour = 'idade_cat=not so young'),
            size = 1)+
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C2, colour = 'idade_cat=old'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C3, colour = 'idade_cat=young'),
            size = 1)


#TRAT
dados_plot = data.frame(x = dados$tempo/365,
                        C1 = logit_mu(optim_result$par[2]) + 
                          (1 - logit_mu(optim_result$par[2]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7])),
                        C2 = logit_mu(optim_result$par[2] + optim_result$par[3]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[3]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[8])),
                        C3 = logit_mu(optim_result$par[2] + optim_result$par[4]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[4]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[9])),
                        C4 = logit_mu(optim_result$par[2] + optim_result$par[5]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[5]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[10])),
                        C5 = logit_mu(optim_result$par[2] + optim_result$par[6]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[6]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[11]))
)

jskm(survfit(Surv(tempo/365, cens) ~ trat, data = dados_tg)) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C1, colour = 'trat=Info ausente'),
            size = 1)+
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C2, colour = 'trat=Medicamento'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C3, colour = 'trat=Medicamento e Radioterapia'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C4, colour = 'trat=Nenhum'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C5, colour = 'trat=Radioterapia'),
            size = 1)

#estagio_p
dados_plot = data.frame(x = dados$tempo/365,
                        C1 = logit_mu(optim_result$par[2]) + 
                          (1 - logit_mu(optim_result$par[2]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7])),
                        C2 = logit_mu(optim_result$par[2] + optim_result$par[3]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[3]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[8])),
                        C3 = logit_mu(optim_result$par[2] + optim_result$par[4]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[4]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[9])),
                        C4 = logit_mu(optim_result$par[2] + optim_result$par[5]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[5]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[10])),
                        C5 = logit_mu(optim_result$par[2] + optim_result$par[6]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[6]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[11]))
)

jskm(survfit(Surv(tempo/365, cens) ~ estagio_p, data = dados_tg)) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C1, colour = 'estagio_p=S1'),
            size = 1)+
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C2, colour = 'estagio_p=S2'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C3, colour = 'estagio_p=S3'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C4, colour = 'estagio_p=S4'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C5, colour = 'estagio_p=SX'),
            size = 1)

#estagio_t
dados_plot = data.frame(x = dados$tempo/365,
                        C1 = logit_mu(optim_result$par[2]) + 
                          (1 - logit_mu(optim_result$par[2]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[6])),
                        C2 = logit_mu(optim_result$par[2] + optim_result$par[3]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[3]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[6] + optim_result$par[7])),
                        C3 = logit_mu(optim_result$par[2] + optim_result$par[4]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[4]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[6] + optim_result$par[8])),
                        C4 = logit_mu(optim_result$par[2] + optim_result$par[5]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[5]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[6] + optim_result$par[9]))
)

jskm(survfit(Surv(tempo/365, cens) ~ estagio_t, data = dados_tg)) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C1, colour = 'estagio_t=T1'),
            size = 1)+
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C2, colour = 'estagio_t=T2'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C3, colour = 'estagio_t=T3'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C4, colour = 'estagio_t=T4'),
            size = 1)

#estagio_n
dados_plot = data.frame(x = dados$tempo/365,
                        C1 = logit_mu(optim_result$par[2]) + 
                          (1 - logit_mu(optim_result$par[2]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7])),
                        C2 = logit_mu(optim_result$par[2] + optim_result$par[3]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[3]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[8])),
                        C3 = logit_mu(optim_result$par[2] + optim_result$par[4]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[4]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[9])),
                        C4 = logit_mu(optim_result$par[2] + optim_result$par[5]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[5]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[10])),
                        C5 = logit_mu(optim_result$par[2] + optim_result$par[6]) + 
                          (1 - logit_mu(optim_result$par[2] + optim_result$par[6]))*
                          sgompertz(dados$tempo/365, optim_result$par[1],
                                    exp(optim_result$par[7] + optim_result$par[11]))
)

jskm(survfit(Surv(tempo/365, cens) ~ estagio_n, data = dados_tg)) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C1, colour = 'estagio_n=N0'),
            size = 1)+
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C2, colour = 'estagio_n=N1'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C3, colour = 'estagio_n=N2'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C4, colour = 'estagio_n=N3'),
            size = 1) +
  geom_line(data=dados_plot, 
            mapping=aes(x=x, y=C5, colour = 'estagio_n=NX'),
            size = 1)

