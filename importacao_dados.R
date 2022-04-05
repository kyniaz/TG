# library('dplyr')
# library('survival')
# library('survminer')
# #library('flexsurvcure')
# 
# clinical <- read.delim("C:/Users/oandr/Downloads/clinical.project-TCGA-BRCA.2021-12-08.tar/clinical.project-TCGA-BRCA.2021-12-08/clinical.tsv")
# 
# dados_surv = clinical |> select(case_id, age_at_diagnosis, gender, race, vital_status, days_to_death, days_to_last_follow_up,
#                                  icd_10_code, treatment_or_therapy, treatment_type, ajcc_pathologic_m, ajcc_pathologic_n, ajcc_pathologic_t, ajcc_pathologic_stage)
# 
# #tratando days to death
# 
# dados_surv$tempo = as.numeric(ifelse(dados_surv$days_to_death == "'--", 
#                                      dados_surv$days_to_last_follow_up, dados_surv$days_to_death))
# 
# remove_aspas = function(i){
#   return(ifelse(dados_surv[,i] == "'--", NA, dados_surv[,i]))
# }
# 
# dados_filtrados = lapply(1:ncol(dados_surv), remove_aspas)
# 
# dados_filtrados = data.frame(dados_filtrados)
# colnames(dados_filtrados) = colnames(dados_surv)
# 
# dados_filtrados$age_at_diagnosis = as.numeric(dados_surv$age_at_diagnosis)
# 
# #convertendo a censura;
# dados_filtrados$cens = ifelse(dados_filtrados$vital_status == 'Alive',0,1)
# 
# dados_filtrados <- dados_filtrados |>
#   mutate_if(sapply(dados_filtrados, is.character), as.factor)
# 
# dados_filtrados_final = dados_filtrados[complete.cases(dados_filtrados |> select(-days_to_last_follow_up, -days_to_death)), ]
# 
# dados_filtrados_final = dados_filtrados_final[dados_filtrados_final$tempo > 0,]
# 
# casos_unicos = unique(dados_filtrados_final$case_id)
# 
# dados_filtrados_final$trat = character(nrow(dados_filtrados_final))
# 
# for (i in 1:length(casos_unicos)){
#   tratamentos = dados_filtrados_final[dados_filtrados_final$case_id == casos_unicos[i],] |>
#     select(treatment_type, treatment_or_therapy)
#   
#   trat = ''
#   if(tratamentos$treatment_or_therapy[1] == 'yes') trat = paste0(trat, tratamentos$treatment_type[1], ' ')
#   if(tratamentos$treatment_or_therapy[2] == 'yes') trat = paste0(trat, tratamentos$treatment_type[2], ' ')
#   
#   if(tratamentos$treatment_or_therapy[1] == 'not reported') trat = paste0(trat, 'Sem Info')
#   if(tratamentos$treatment_or_therapy[2] == 'not reported') trat = paste0(trat, 'Sem Info')
#   
#   dados_filtrados_final[dados_filtrados_final$case_id == casos_unicos[i],]$trat = trat
# }
# 
# dados_filtrados_final$trat = ifelse(dados_filtrados_final$trat == '', 'Nenhum', dados_filtrados_final$trat)
# dados_filtrados_final$trat = ifelse(dados_filtrados_final$trat %in% c('Pharmaceutical Therapy, NOS Radiation Therapy, NOS ',
#                                     'Radiation Therapy, NOS Pharmaceutical Therapy, NOS ')
#                                     , 'Medicamento e Radioterapia', dados_filtrados_final$trat)
# dados_filtrados_final$trat = ifelse(dados_filtrados_final$trat %in% c('Pharmaceutical Therapy, NOS ',
#                                                                       'Pharmaceutical Therapy, NOS Sem Info')
#                                     , 'Medicamento', dados_filtrados_final$trat)
# dados_filtrados_final$trat = ifelse(dados_filtrados_final$trat %in% c('Radiation Therapy, NOS ',
#                                                                       'Radiation Therapy, NOS Sem Info')
#                                     , 'Radioterapia', dados_filtrados_final$trat)
# dados_filtrados_final$trat = ifelse(dados_filtrados_final$trat %in% c('Sem InfoSem Info',
#                                                                       'Sem Info')
#                                     , 'Info ausente', dados_filtrados_final$trat)
# table(dados_filtrados_final$trat)
# 
# dados_filtrados_final = dados_filtrados_final |> select(-treatment_type, - treatment_or_therapy) |> unique()
# 
# dados_filtrados_final$estagio_m = as.factor(ifelse(dados_filtrados_final$ajcc_pathologic_m %in% c('cM0 (i+)'), 
#                                                    "M0", as.character(dados_filtrados_final$ajcc_pathologic_m)))
# 
# dados_filtrados_final$estagio_n = ifelse(dados_filtrados_final$ajcc_pathologic_n %in% c('N0', 'N0 (i-)','N0 (i+)','N0 (mol+)'), 'N0',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_n %in% c('N1', 'N1a', 'N1b', 'N1c','N1mi'),'N1',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_n %in% c('N2', 'N2a'),'N2',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_n %in% c('N3', 'N3a','N3c'),'N3',
#                                          'NX'))))
# 
# dados_filtrados_final$estagio_t = ifelse(dados_filtrados_final$ajcc_pathologic_t %in% c('T1', 'T1a','T1b','T1c'), 'T1',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_t %in% c('T2', 'T2a', 'T2b'),'T2',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_t %in% c('T3', 'T3a'),'T3',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_t %in% c('T4', 'T4b','T4d'),'T4',
#                                         'T4'))))
# 
# table(dados_filtrados_final$ajcc_pathologic_stage)
# 
# dados_filtrados_final$estagio_p = ifelse(dados_filtrados_final$ajcc_pathologic_stage %in% c('Stage I', 'Stage IA','Stage IB'), 'S1',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_stage %in% c('Stage II', 'Stage IIA', 'Stage IIB'),'S2',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_stage %in% c('Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC'),'S3',
#                                          ifelse(dados_filtrados_final$ajcc_pathologic_stage %in% c('Stage IV'),'S4',
#                                         'SX'))))
# table(dados_filtrados_final$estagio_p)
# 
# dados_tg = dados_filtrados_final |> select (case_id, age_at_diagnosis, gender, race, vital_status, days_to_death,
#                                             days_to_last_follow_up, icd_10_code, ajcc_pathologic_m, ajcc_pathologic_n,
#                                             ajcc_pathologic_t, ajcc_pathologic_stage, tempo, cens, trat, estagio_m, estagio_n, estagio_p, estagio_t)
# 
# write.csv(dados_tg, file = 'dados_tg.csv')
