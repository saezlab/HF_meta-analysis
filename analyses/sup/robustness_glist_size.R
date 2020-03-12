# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: This script that shows the robustness of using different
#' gene list sizes to calculate disease scores and enrichment scores

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

experiments = names(METAheart)
names(experiments) = experiments

t_matrix = get_all_limma(meta_list = METAheart,
                         limma_column = "t")

# 1. Robustness of disease score

n_genes = c(50,100,200,500,1000)
names(n_genes) = n_genes

pairwise_auc = lapply(n_genes, pairwise_ds, 
                      experiments = experiments,
                      meta_list = METAheart,
                      t_matrix = t_matrix) 

pairwise_auc = enframe(pairwise_auc,name = "n_genes") %>% unnest()

saveRDS(pairwise_auc, file =  "data/figure_objects/robust_ds.rds")

# 2. Robustness of ES

pairwise_es_res = lapply(n_genes, pairwise_ES,
                      meta_list = METAheart) 

pairwise_es_res = enframe(pairwise_es_res,
                          name = "n_genes") %>% unnest()

saveRDS(pairwise_es_res, file =  "data/figure_objects/robust_es.rds")

