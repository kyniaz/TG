library('dplyr')
library('survival')
library('survminer')
library('lubridate')
#library('flexsurvcure')
library('ggfortify')

library('Rcpp')
sourceCpp(file = 'gompertz.cpp')


clinical <- read.delim("C:/Users/oandr/Downloads/clinical.project-TCGA-BRCA.2021-12-08.tar/clinical.project-TCGA-BRCA.2021-12-08/clinical.tsv")

#clinical <- read.delim("C:/Users/oandr/Downloads/clinical.project-TCGA-LUAD.2021-12-08.tar/clinical.tsv", header=T)

dados_surv = clinical %>% select(case_id, age_at_diagnosis, gender, race, vital_status, days_to_death, days_to_last_follow_up,
                                 icd_10_code, treatment_or_therapy, treatment_type, ajcc_pathologic_m, ajcc_pathologic_n, ajcc_pathologic_t, ajcc_pathologic_stage)

#tratando days to death

dados_surv$tempo = as.numeric(ifelse(dados_surv$days_to_death == "'--", 
                                     dados_surv$days_to_last_follow_up, dados_surv$days_to_death))

remove_aspas = function(i){
  return(ifelse(dados_surv[,i] == "'--", NA, dados_surv[,i]))
}

dados_filtrados = lapply(1:ncol(dados_surv), remove_aspas)

dados_filtrados = data.frame(dados_filtrados)
colnames(dados_filtrados) = colnames(dados_surv)

dados_filtrados$age_at_diagnosis = as.numeric(dados_surv$age_at_diagnosis)

#convertendo a censura;
dados_filtrados$cens = ifelse(dados_filtrados$vital_status == 'Alive',0,1)

dados_filtrados <- dados_filtrados %>%
  mutate_if(sapply(dados_filtrados, is.character), as.factor)

#combinar tratamentos...

dados_filtrados_final = dados_filtrados[complete.cases(dados_filtrados %>% select(-days_to_last_follow_up, -days_to_death)), ]

dados_filtrados_final = dados_filtrados_final[dados_filtrados_final$tempo > 0,]

dados_filtrados_final = dados_filtrados_final %>% select(-treatment_type, - treatment_or_therapy) %>% unique()

dados_filtrados_final = dados_filtrados_final[!(dados_filtrados_final$race == 'american indian or alaska native'),]

######################
#TTT plot com censura#
######################

tempo = dados_filtrados_final$tempo
censura = dados_filtrados_final$cens

o=order(tempo)
t=tempo[o]
cens=censura[o]


n=length(t)
r=sum(cens)

j=1
TF=numeric()
MON=numeric()
I=numeric()
TTT=numeric()
Fi=numeric()
F_var=numeric()
S=numeric()

TF[j]=0
MON[j]=0
F_var[j]=0
S[j]=1
TTT[j]=0
i=1

while(i<(n+1)){
  if(cens[i]==1){
    j=j+1
    TF[j]=t[i]
    NI=n-i+1
    I=((n+1)-MON[j-1])/(1+NI)
    MON[j]=MON[j-1]+I
    F_var[j]=MON[j]/n
    S[j]=1-F_var[j]
  }
  i=i+1
}

TF[r+2]=t[n]
F_var[r+2]=1
TTT[1]=0

for(j in 2:(r+2)){
  TTT[j]=TTT[j-1]+n*S[j-1]*(TF[j]-TF[j-1])
}

for(j in 1:(r+2)){
  Fi[j]=TTT[j]/TTT[r+2]
}

ggplot(data = data.frame(F_var, Fi), aes(x = F_var, y = Fi))+
  geom_point() +
  theme_minimal() + 
  labs(x = "F(t)", y = "TTT com censuras") +
  geom_abline(slope=1, intercept=0) +
  geom_line() + 
  xlim(0,1) + ylim(0,1)


kaplan_meier = survfit(Surv(tempo/365, cens) ~ 1, data = dados_filtrados_final)

plot_dados = data.frame(kaplan_meier$time, kaplan_meier$surv, kaplan_meier$n.event)
colnames(plot_dados) = c('Tempo', 'Sobrevivência', 'Evento')

ggplot() +
  geom_line(aes(x = Tempo, y = Sobrevivência), data = plot_dados, size = 1, color = 'steelblue') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #geom_point(aes( x = Tempo, y = Sobrevivência), data = plot_dados[plot_dados$Evento == 1,], shape = 3, color = 'brown2') +
  labs(x = 'Tempo em anos') +
  ylim(0, 1) + 
  theme_minimal() 

