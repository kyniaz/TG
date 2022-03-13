library('dplyr')
library('survival')
library('survminer')
#library('flexsurvcure')

library('Rcpp')
sourceCpp(file = 'gompertz.cpp')


clinical <- read.delim("C:/Users/oandr/Downloads/clinical.project-TCGA-BRCA.2021-12-08.tar/clinical.project-TCGA-BRCA.2021-12-08/clinical.tsv")

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

dados_filtrados_final = dados_filtrados[complete.cases(dados_filtrados %>% select(-days_to_last_follow_up, -days_to_death)), ]

dados_filtrados_final = dados_filtrados_final[dados_filtrados_final$tempo > 0,]

dados_filtrados_final = dados_filtrados_final %>% select(-treatment_type, - treatment_or_therapy) %>% unique()

dados_filtrados_final = dados_filtrados_final[!(dados_filtrados_final$race == 'american indian or alaska native'),]


