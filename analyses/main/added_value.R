# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we test the differences in the predictions between the signature observed in our meta-analysis
#' and individual signatures from experiments

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(PLIER)
library(reshape2)
library(ROCR)
library(fgsea)
library(dplyr)
library(tidyr)
library(WriteXLS)

experiments = names(METAheart)
names(experiments) = experiments

# For labeling
t_matrix = get_all_limma(meta_list = METAheart,
                         limma_column = "t")

# Here we get AUC for pairwise classifiers
pairwise_500 = pairwise_ds(experiments = experiments,
                           meta_list = METAheart,
                           t_matrix = t_matrix,
                           ngenes = 500) #Second page excel

# Here we get the results from performing the meta-analysis

fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10)

sum(fisher_rank < .00005)

genes = names(fisher_rank)

# Here we get the calculations of using the top N genes from the meta-analysis

ds_top = getRisk_Stats_v2(Experiment_List = METAheart,
                          limma_t_mat = t_matrix, 
                          genes = names(fisher_rank[1:500]))

ds_top_predictions = enframe(lapply(ds_top, function(x) {
  enframe(x[["SingleAUC"]])
})) %>% unnest()

colnames(ds_top_predictions) = c("PredictedExperiment", 
                                 "PredictorExperiment", 
                                 "meta_auc")

# Here we merge them 

comp_df = left_join(pairwise_500,
                    ds_top_predictions) %>%
          dplyr::filter(PredictedExperiment != PredictorExperiment) %>%
  dplyr::select(PredictorExperiment, 
                PredictedExperiment,
                single, meta_auc)

print("Are AUCs better?")

wilcox.test(comp_df$meta_auc,comp_df$single,paired = T,alternative = "greater")

comp_df = comp_df %>%
  gather("type","AUROC",
         -PredictorExperiment,
         -PredictedExperiment)

saveRDS(comp_df, "data/figure_objects/added_value_ds.rds")


# Here we test the added value at consistency

study_deg_list_up = lapply(METAheart, function(x){
  deg = dplyr::filter(x$HF_limma,
                     ID %in% names(fisher_rank[1:500])) %>%
    filter(t >0)
  return(deg[[1]])
})

study_deg_list_down = lapply(METAheart, function(x){
  deg = dplyr::filter(x$HF_limma,
                      ID %in% names(fisher_rank[1:500])) %>%
    filter(t < 0)
  return(deg[[1]])
})

# upregulation

up_ES = lapply(experiments, function(x){
  
  stat_rank = METAheart[[x]][["HF_limma"]][["t"]]
  names(stat_rank) = METAheart[[x]][["HF_limma"]][["ID"]]
  stat_rank = sort(stat_rank)
  set.seed(1234)
  
  up_row = as_tibble(fgsea(pathways = study_deg_list_up,
                           stats = stat_rank,nperm = 1000)) %>%
    dplyr::select(pathway,ES)
})

up_ES = up_ES %>% 
  enframe("Reference") %>% unnest()

colnames(up_ES) = c("Reference","DEG","ES")


# downregulation

down_ES = lapply(experiments, function(x){
  
  stat_rank = METAheart[[x]][["HF_limma"]][["t"]]
  names(stat_rank) = METAheart[[x]][["HF_limma"]][["ID"]]
  stat_rank = sort(stat_rank)
  set.seed(1234)
  
  up_row = as_tibble(fgsea(pathways = study_deg_list_down,
                           stats = stat_rank,nperm = 1000)) %>%
    dplyr::select(pathway,ES)
})

down_ES = down_ES %>% 
  enframe("Reference") %>% unnest()

colnames(down_ES) = c("Reference","DEG","ES")

# Previous results
up_ES_single = readRDS(file = "data/figure_objects/up_ES.rds")
down_ES_single = readRDS(file = "data/figure_objects/down_ES.rds")


# Merge them
print("Is the ES better?")

comp_ES_up= left_join(up_ES_single,
                      up_ES, by = c("Reference","DEG")) %>%
  dplyr::filter(Reference != DEG) 

wilcox.test(comp_ES_up$ES.y,
            comp_ES_up$ES.x,
            paired = T,
            alternative = "greater")

comp_ES_down = left_join(down_ES_single,
                      down_ES, by = c("Reference","DEG")) %>%
  dplyr::filter(Reference != DEG)

wilcox.test(comp_ES_down$ES.y,
            comp_ES_down$ES.x,
            paired = T,
            alternative = "less")

comp_ES = bind_rows("upregulated"=comp_ES_up, 
                    "downregulated" =comp_ES_down,
                    .id = "comparison") %>% 
  mutate(single = ES.x, 
         meta_analysis = ES.y) %>%
  dplyr::select(c("comparison","Reference","DEG",
                  "single","meta_analysis")) %>%
  gather("type","ES", 
         -comparison,-Reference, 
         -DEG)

saveRDS(comp_ES, "data/figure_objects/added_value_es.rds")

  












