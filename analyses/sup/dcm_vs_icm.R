# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: This script estimates the variability accounted for
#' diverse etiologies of heart failure

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")

library(PLIER)
library(sjstats)
library(tidyverse)

METAheart = readRDS(file = "data/METAheart.rds") #main object
covariate_summary =readRDS(file = "data/covariate_summary.rds")

#0. Preparing data

# Gene centering: we expect that ICM and DCM interact in the same direction...
# It may overestimate differences, but good in practice

dcm_studies = covariate_summary %>% 
              dplyr::filter(DCM == "yes") %>%
              dplyr::select("ID")

dcm_studies = dcm_studies[[1]]

load("data/dictionaryIDs.ro")
new_ids = dictionary %>% dplyr::filter(GEO_ID %in% dcm_studies)
dcm_studies = new_ids$newID

disease_METAheart = lapply(METAheart[dcm_studies], function(x){
  
  x$TARGETS = dplyr::filter(x$TARGETS,
                HeartFailure == "yes")
  
  x$GEX = x$GEX[,x$TARGETS$Sample]
  
  x[["GEX_norm"]] = PLIER::rowNorm(x[["GEX"]])
  
  return(x)
  
})

disease_targets = get_tibble_union(disease_METAheart,"TARGETS") %>% 
  dplyr::select(Sample,ExpID,DCM) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_"))

meta_gex_scale = get_complete_gex(meta_list = disease_METAheart,
                                  complete_targets = disease_targets,
                                  gex_key = "GEX_norm")

#1. Performing a PCA

meta_gex_naomit_scale = meta_gex_scale[rowSums(is.na(meta_gex_scale))==0,]
pca_meta_scale = prcomp(t(meta_gex_naomit_scale[,disease_targets$grl_id]),
                        center = T,scale. = T)
pca_meta_sum_scale = summary(pca_meta_scale)

pca_meta_sum_scale[["plotdf"]] = disease_targets %>% 
                                mutate("PC1" = pca_meta_scale$x[,1],
                                       "PC2" = pca_meta_scale$x[,2])

saveRDS(pca_meta_sum_scale, 
        file = "data/figure_objects/dcm_icm_pca.rds")


#2. Fitting a linear model to each PC to find associations with ICM or DCM

pcs_study = run_anovastats_single(numeric_matrix = t(pca_meta_sum_scale$x),
                                  targets = disease_targets,
                                  factor_a = "DCM",
                                  pval = 0.05)

pcs_study = mutate(pcs_study,
                   prop_var = pca_meta_sum_scale$importance[2,pcs_study$PC]) %>%
            arrange(desc(prop_var)) %>% 
            select(PC, factor_a_vect, prop_var) %>% 
            mutate(factor_a_vect = as.character(factor_a_vect))

# ugly trick to make table
total_row = c("TOTAL","",sum(pcs_study$prop_var))
pcs_study = pcs_study %>% mutate(prop_var = as.character(prop_var))
names(total_row) = colnames(pcs_study) = c("PCs","p_value","Prop. Var")


pcs_study = bind_rows(pcs_study, total_row)

print(pcs_study,n=100)

saveRDS(pcs_study, file = "data/figure_objects/dcm_icm_pcs.rds")












