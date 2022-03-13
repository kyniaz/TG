############# 
############ Gompertz sem cura, com Cura.

#Com mistura e sem mistura, bayesiano apenas..

source('importacao_dados.R')
library('microbenchmark')
library('JuliaCall')

julia_setup()
####Sem cura -----

par(mfrow = c(1,1))

plot(survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final), 
     xlab = "Anos", 
     ylab = "Sobrevivência")

metropolis_gibbs <- function(delta, a, b, rprop, ldprop, tempos, d, B = 10^3, start = 0.5)
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
      
      if(is.nan(lratio)) chain[[ii]] <- chain[[ii-1]]
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

gamma_hiper_a = c(1,1)
gamma_hiper_b = c(4,1)

log_veros = function(tempos, cens, alfa, beta_p){
  return(sum(cens*dgompertz(tempos, a = alfa, b = beta_p, log_opt = T) + 
               (1-cens)*(sgompertz(tempos, a = alfa, b = beta_p,log_opt = T))))
}


ldtgt1 <- function(tempos, d, a, b) {
  #if(prob>1) {print(paste0('delta(1) :', delta, ' e prob: ', prob))}
  return(log_veros(tempos, d, a, b) + dgamma(a,gamma_hiper_a[1],gamma_hiper_a[2], log = T) + dgamma(b,gamma_hiper_b[1],gamma_hiper_b[2], log = T))
}

ldtgt2 <- function(tempos, d, a, b) {
  return(log_veros(tempos, d, a, b) + dgamma(a,gamma_hiper_a[1],gamma_hiper_a[2], log = T) + dgamma(b,gamma_hiper_b[1],gamma_hiper_b[2], log = T))
}

####### Propostas ----

rprop_gamma = function(ant) rgamma(1, ant, 1)
ldprop_gamma = function(ant, prop) dgamma(prop, ant, 1, log = T)
#R

metropolis_R = function(tempos, eventos, B , BB, start_a, start_b){
  a = rep(NA, BB)
  b = rep(NA, BB)
  
  a[1] = 0.1
  b[1] = 3.5
  
  tempos = dados_filtrados_final$tempo/365
  eventos = dados_filtrados_final$cens
  
  for (ii in 2:BB){
    delta = trunc(runif(1, 1, 3))
    a[ii] = ifelse(delta == 1, metropolis_gibbs(delta, a[ii-1], b[ii-1], rprop_gamma, 
                                                ldprop_gamma, tempos, eventos , B, start = start_a), a[ii-1])
    b[ii] = ifelse(delta == 2, metropolis_gibbs(delta, a[ii-1], b[ii-1], rprop_gamma, 
                                                ldprop_gamma, tempos, eventos , B, start = start_b), b[ii-1])
  }
  return(list(a = a, b = b))
}


strt = Sys.time()
metropolis_R(tempos, eventos, 500, 200, 0.05, 4)
Sys.time() - strt

julia_source('non_cure.jl')

###### Com Cura - Mistura;

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


############### C
