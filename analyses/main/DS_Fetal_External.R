# Author: Jan Lanzer 2019
# Description: This script performs the diseaes score (risk score) calculation for fetal and 
library(tidyverse)
library(WriteXLS)

# load sorurce   functions
source("src/data_utils.R") 
source("src/misc_utils.R")

# load objects containing data from heart failure studies
METAheart= readRDS(file = "data/METAheart.rds")
external_experiments = readRDS("data/external_METAheart.rds")
fetal_experiments = readRDS("data/fetal_METAheart.rds")
meta_rank = readRDS("data/shiny/fisher_rank.rds")



#### Part 1. Analysis of External Datasets ####

# A) Disease Score calculation
# calculate coefficients (t-values)
coef_matrix = get_all_limma(meta_list = METAheart,limma_column = "t")

# Predicting genes (top 500)
predictor_genes = names(meta_rank[1:500])

# calculating the disease score for external_experiments
mat_ds = getRisk_Stats_v2(Experiment_List = external_experiments,
                           genes = predictor_genes,
                           limma_t_mat = coef_matrix)

# calculate the mean disease score for each external_experiment 
mat_ds_means = lapply(mat_ds , function(x){
  data.frame("Risk_Score" = rowMeans(x$RiskMatrix))
    })

# add sample information from target file to disease scores
for (study in names(external_experiments)){
  mat_ds_means[[study]] = mat_ds_means[[study]]%>%
    rownames_to_column("Sample") %>% 
    as_tibble() %>% 
    full_join(external_experiments[[study]]$TARGETS) %>%
    select(Sample, Risk_Score, HeartFailure) %>% 
    mutate(HeartFailure= as.factor(HeartFailure))
}

# save disease score
saveRDS(mat_ds_means,
        file = "data/paper_sup/ds_external_experiment.rds")
saveRDS(mat_ds,
        file = "data/paper_sup/ds_external_experiment_AUC.rds")


# B) Disease Score modelling
# logistic regression (LR) - to be deleted, not a good model here.
# create empty data frame to store LR results
nstudies =  length(names(mat_ds_means))
logistic_results =data.frame("Coef"= rep(NA, nstudies),"Pval"= rep(NA,nstudies))
rownames(logistic_results) = names(mat_ds_means)

# perform LR
for (study in rownames(logistic_results)){
  glm.fit = glm(HeartFailure ~ Risk_Score, data = mat_ds_means[[study]], family  =  binomial)
  logistic_results[study,"Pval"]= summary(glm.fit)$coef["Risk_Score","Pr(>|z|)"]
  logistic_results[study,"Coef"]= summary(glm.fit)$coef["Risk_Score","Estimate"]
  }

# linear model (LM) 
# create empty data frame to store LM results
nstudies =  length(names(mat_ds_means))
lm_results =data.frame("Coef"= rep(NA, nstudies),"Pval"= rep(NA,nstudies))
rownames(lm_results) = names(mat_ds_means)

#perform LM
for (study in names(mat_ds_means)){
  lm.ds= lm(Risk_Score ~ HeartFailure, data= mat_ds_means[[study]])
  lm_results[study,"Pval"]= summary(lm.ds)$coef["HeartFailureyes","Pr(>|t|)"]
  lm_results[study,"Coef"]= summary(lm.ds)$coef["HeartFailureyes","Estimate"]
}

lm_results = lm_results %>% rownames_to_column("Study")

# save results of LM
saveRDS(lm_results,
        file = "data/paper_sup/ds_external_experiment_lm.rds")

# save results of DS as .xlsx file
ds_GSE84796  = mat_ds_means$GSE84796
ds_GSE4172   = mat_ds_means$GSE4172
ds_GSE9800   = mat_ds_means$GSE9800
ds_GSE10161  = mat_ds_means$GSE10161
ds_GSE3586   = mat_ds_means$GSE3586

WriteXLS(x = c("ds_GSE84796","ds_GSE4172","ds_GSE9800","ds_GSE10161","ds_GSE3586", "lm_results"),
         ExcelFileName = "data/paper_sup/External_DS.xlsx",
         SheetNames = c("ds_GSE84796","ds_GSE4172","ds_GSE9800","ds_GSE10161","ds_GSE3586","lm_results")
)


#### Part 2. Analysis of fetal datasets ########################################

# renaming trick (rename fetal to HeartFailure) to run functions from source that requires HeartFailure category
fetal_experiments_copy= fetal_experiments
  for (study in names(fetal_experiments_copy)){
    fetal_experiments_copy[[study]]$TARGETS = fetal_experiments_copy[[study]]$TARGETS %>%
      rename(HeartFailure_original = HeartFailure) %>%
      rename(HeartFailure = fetal)
  }


# A) Disease Score calculation
# calculating the disease score for fetal experiments
mat_ds = getRisk_Stats_v2(Experiment_List = fetal_experiments_copy,
                          genes = predictor_genes,
                          limma_t_mat = coef_matrix)

# calculate the mean disease score for each fetal experiment 
mat_ds_means = lapply(mat_ds , function(x){
  data.frame("Risk_Score" = rowMeans(x$RiskMatrix))
})

# add sample information from target file to disease scores
for (study in names(fetal_experiments_copy)){
  mat_ds_means[[study]] = mat_ds_means[[study]]%>%
    rownames_to_column("Sample") %>% 
    as_tibble() %>% 
    full_join(fetal_experiments_copy[[study]]$TARGETS) %>%
    select(Sample, Risk_Score, HeartFailure) %>% 
    mutate(HeartFailure= as.factor(HeartFailure))
}

#Repeat disease score calculation for PRJNA522417
#Reason: For the fetal study spurell19, we cannot use the t-value matrix containing the very same study to predict itself, this 
#would be overfitting and therefore they have to be taken out before calculating disease score

# A) Disease Score calculation
# calculating the disease score for fetal experiments
mat_ds_spurrell = getRisk_Stats_v2(Experiment_List = fetal_experiments_copy,
                                   genes = predictor_genes,
                                   limma_t_mat = coef_matrix[,colnames(coef_matrix) != "PRJNA522417"]) #  here we remove the study from the coef matrix

# calculate the mean disease score for each fetal experiment 
mat_ds_means_spurrell = lapply(mat_ds_spurrell , function(x){
  data.frame("Risk_Score" = rowMeans(x$RiskMatrix))
})

# add sample information from target file to disease scores
for (study in names(fetal_experiments_copy)){
  mat_ds_means_spurrell[[study]] = mat_ds_means_spurrell[[study]]%>%
    rownames_to_column("Sample") %>% 
    as_tibble() %>% 
    full_join(fetal_experiments_copy[[study]]$TARGETS) %>%
    select(Sample, Risk_Score, HeartFailure) %>% 
    mutate(HeartFailure= as.factor(HeartFailure))
}


#overwriting the overfitted ds with the correct values
mat_ds_means$PRJNA522417 = mat_ds_means_spurrell$PRJNA522417
mat_ds$PRJNA522417 = mat_ds_spurrell$PRJNA522417


# save disease score
saveRDS(mat_ds_means,
        file = "data/paper_sup/ds_fetal_experiment.rds")
saveRDS(mat_ds,
        file = "data/paper_sup/ds_fetal_experiment_AUC.rds")


# B) Disease Score modelling
# linear model (LM) 
# create empty data frame to store LM results
nstudies =  length(names(mat_ds_means))
lm_results =data.frame("Coef"= rep(NA, nstudies),"Pval"= rep(NA,nstudies))
rownames(lm_results) = names(mat_ds_means)

#perform LM
for (study in names(mat_ds_means)){
  lm.ds= lm(Risk_Score ~ HeartFailure, data= mat_ds_means[[study]])
  lm_results[study,"Pval"]= summary(lm.ds)$coef["HeartFailureyes","Pr(>|t|)"]
  lm_results[study,"Coef"]= summary(lm.ds)$coef["HeartFailureyes","Estimate"]
}

lm_results = lm_results %>% rownames_to_column("Study")

# save results of LM
saveRDS(lm_results,
        file = "data/paper_sup/ds_fetal_experiment_lm.rds")

# save results of DS as .xlsx file
ds_GSE52501 = mat_ds_means$GSE52601
ds_PRNJA522417 = mat_ds_means$PRJNA522417
WriteXLS(x= c("ds_GSE52501","ds_PRNJA522417","lm_results"),
         ExcelFileName = "data/paper_sup/Fetal_DS.xlsx",
         SheetNames = c("ds_GSE52501","ds_PRNJA522417","lm_results")
)








