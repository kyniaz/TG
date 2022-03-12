par(mfrow = c(2,2))

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

############### GOMPERTZ ----

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dgompertz(t, a = par[1], b = par[2], log_opt = T) + (1-d)*log(par[3] + (1-par[3])*(1-pgompertz(t, a = par[1], b = par[2])))))
}

param = c(1,1,0.5)


teste = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste$par

curve(teste$par[3] + (1-teste$par[3])*(1-pgompertz(x, teste$par[1], teste$par[2])), col = 'red', add = T)


log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dgompertz(t, a = par[1], b = par[2], log_opt = T) + log(par[3])*pgompertz(t, a = par[1], b = par[2]) ))
}

param2 = c(1,1,0.3)

teste2 = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste2$par 


curve((teste2$par[3])^(pgompertz(x, a = teste2$par[1], b = teste2$par[2])), col = 'blue', add = T)

############### GAMMA -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dgamma(t, par[1], par[2], log = T) + 
               (1-d)*log(par[3] + (1-par[3])*(1-pgamma(t, par[1], par[2])))))
}

param = c(1,1,0.5)


teste = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste$par

curve(teste$par[3] + (1-teste$par[3])*(1-pgamma(x, teste$par[1], teste$par[2])), col = 'red', add = T)


log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dgamma(t, par[1], par[2], log = T) + 
               log(par[3])*pgamma(t, par[1], par[2]) ))
}

param2 = c(1,1,0.3)

teste2 = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste2$par 


curve((teste2$par[3])^(pgamma(x, teste2$par[1], teste2$par[2])), col = 'blue', add = T)

########## Exponencial -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[2])) + d*dexp(t, par[1], log = T) + 
               (1-d)*log(par[2] + (1-par[2])*(1-pexp(t, par[1])))))
}

param = c(0.01,0.5)

teste = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste$par

curve(teste$par[2] + (1-teste$par[2])*(1-pexp(x, teste$par[1])), col = 'red', add = T)

log_veros2 = function(par){
  return(sum(d*log(-log(par[2])) + d*dexp(t, par[1], log = T) + 
               log(par[2])*pexp(t,par[1]) ))
}

param2 = c(0.01,0.3)

teste2 = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste2$par 

curve((teste2$par[2])^(pexp(x, teste2$par[1])), col = 'blue', add = T)

########## WEIBULL -----

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

t = dados_filtrados_final$tempo/365
d = dados_filtrados_final$cens

log_veros = function(par){
  return(sum(d*(log(1-par[3])) + d*dweibull(t, par[1], par[2], log = T) + 
               (1-d)*log(par[3] + (1-par[3])*(1-pweibull(t, par[1], par[2])))))
}

param = c(1,1,0.5)


teste = optim(param, log_veros, control = list(fnscale = -1, maxit = 500),
              method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste$par

curve(teste$par[3] + (1-teste$par[3])*(1-pweibull(x, teste$par[1], teste$par[2])), col = 'red', add = T)


log_veros2 = function(par){
  return(sum(d*log(-log(par[3])) + d*dweibull(t, par[1], par[2], log = T) + 
               log(par[3])*pweibull(t, par[1], par[2]) ))
}

param2 = c(1,1,0.3)

teste2 = optim(param2, log_veros2, control = list(fnscale = -1, maxit = 500),
               method="L-BFGS-B", lower = c(0.0001,0.0001,0.0001), upper = c(Inf,Inf,0.9999))

teste2$par 


curve((teste2$par[3])^(pweibull(x, teste2$par[1], teste2$par[2])), col = 'blue', add = T)


############# Bayesiana ----

###BAYES

metropolis_gibbs <- function(prob, a, b, rprop, ldprop, t, d, B = 10^4, start = 0.5)
{
  chain <- as.list(rep(NA, B))
  chain[[1]] <- start
  
  if(delta == 1) {
    for(ii in 2:B)
    {
      prop <- rprop(chain[[ii-1]])
      
      lratio <- ldtgt1(t, d, prob, prop, b) -ldtgt1(t, d, prob, chain[[ii-1]], b) +
        ldprop(prop,chain[[ii-1]])-
        ldprop(chain[[ii-1]],prop)
      
      if(is.nan(lratio)) chain[[ii]] = chain[[ii-1]]
      else if(log(runif(1)) <= lratio) chain[[ii]] <- prop
      else{
        chain[[ii]] <- chain[[ii-1]]
      }
    }
  }
  if(delta == 2) {
    for(ii in 2:B)
    {
      prop <- rprop(chain[[ii-1]])
      
      lratio <- ldtgt2(t, d, prob, a, prop) - ldtgt2(t, d, prob, a, chain[[ii-1]]) +
        ldprop(prop,chain[[ii-1]])-
        ldprop(chain[[ii-1]],prop)
      
      if(is.nan(lratio)) chain[[ii]] = chain[[ii-1]]
      else if(log(runif(1)) <= lratio) chain[[ii]] <- prop
      else{
        chain[[ii]] <- chain[[ii-1]]
      }
    }
  }
  if(delta == 3) {
    for(ii in 2:B)
    {
      prop <- rprop(chain[[ii-1]])
      
      lratio <- ldtgt3(t, d, prop, a, b) - ldtgt3(t, d, chain[[ii-1]], a, b) +
        ldprop(prop,chain[[ii-1]])-
        ldprop(chain[[ii-1]],prop)
      
      if(is.nan(lratio)) chain[[ii]] = chain[[ii-1]]
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

gamma_hiper = c(1,1)

log_veros_mix = function(tempos, cens, mix, alfa, beta_p){
  return(sum(cens*(log(1-mix)) + cens*dgompertz(tempos, a = alfa, b = beta_p, log_opt = T) + (1-cens)*log(mix + (1-mix)*(1-pgompertz(tempos, a = alfa, b = beta_p)))))
}

log_veros_nomix = function(tempos, cens, mix, alfa, beta_p){
  return(sum(cens*(log(-log(mix))) + cens*dgompertz(tempos, a = alfa, b = beta_p, log_opt = T) + log(mix)*pgompertz(tempos, a = alfa, b = beta_p)))
}


ldtgt1 <- function(t, d, prob, a, b) {
  #if(prob>1) {print(paste0('delta(1) :', delta, ' e prob: ', prob))}
  return(log_veros_nomix(t, d, prob, a, b) + dgamma(a,gamma_hiper[1],gamma_hiper[2], log = T) + dgamma(b,gamma_hiper[1],gamma_hiper[2], log = T))
}

ldtgt2 <- function(t, d, prob, a, b) {
  return(log_veros_nomix(t, d, prob, a, b) + dgamma(a,gamma_hiper[1],gamma_hiper[2], log = T) + dgamma(b,gamma_hiper[1],gamma_hiper[2], log = T))
}

ldtgt3 <- function(t, d, prob, a, b) {
  return(log_veros_nomix(t, d, prob, a, b) + dgamma(a,gamma_hiper[1],gamma_hiper[2], log = T) + dgamma(b,gamma_hiper[1],gamma_hiper[2], log = T) + log(1))
}
####### Propostas ###################

rprop_beta <- function(ant) rbeta(1, 1, 1)
ldprop_beta <- function(ant,prop) dbeta(prop, 1,1, log = T)

rprop_gamma = function(ant) rexp(1, 1/ant)
ldprop_gamma = function(ant, prop) dexp(prop, 1/ant, log = T)

B = 1000
BB = 1000
a = rep(NA, BB)
b = rep(NA, BB)
prob = rep(NA, BB)

a[1] = 0.5
b[1] = 0.5
prob[1] = 0.5

start_time <- Sys.time()

for (ii in 2:BB){
  delta = trunc(runif(1, 1, 4))
  a[ii] = ifelse(delta == 1, metropolis_gibbs(prob[ii-1], a[ii-1], b[ii-1], rprop_gamma, ldprop_gamma, t, d , B, start = 0.5), a[ii-1])
  b[ii] = ifelse(delta == 2, metropolis_gibbs(prob[ii-1], a[ii-1], b[ii-1], rprop_gamma, ldprop_gamma, t, d , B, start = 0.5), b[ii-1])
  prob[ii] = ifelse(delta == 3, metropolis_gibbs(prob[ii-1], a[ii-1], b[ii-1], rprop_beta, ldprop_beta, t, d , B, start = 0.5), prob[ii-1])
}

end_time <- Sys.time()
end_time - start_time
