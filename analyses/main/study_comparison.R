# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script generates the main analysis of
# reproducibility. Jaccard Index / Disease Score / Enrichment Analysis

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(PLIER)
library(reshape2)
library(ROCR)
library(fgsea)
library(WriteXLS)

# 0. Preparing data
experiments = names(METAheart)
names(experiments) = experiments

# For labeling
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)
# 1. Jaccaard Index

# Extra step use all data first to compare platform coverage
study_genelist_all = lapply(METAheart, function(x){
  gene_list = x$HF_limma[[1]]
  return(gene_list)
}) 

# Make a matrix/data-frame with pairwaise comparisons
jaccard_res_all =  enframe(lapply(experiments, function(set_a){
  genes_a = study_genelist_all[[set_a]]
  j_ix_a = lapply(experiments, function(set_b){
    
    genes_b = study_genelist_all[[set_b]]
    
    #Jaccard Index
    j_ix = length(intersect(genes_a,genes_b))/length(union(genes_a,genes_b))
    
    return(j_ix)
  })
  
  j_ix_a = enframe(j_ix_a, "StudyB","JaccardIx") %>% 
    unnest()
  
  return(j_ix_a)
  
}), "StudyA") %>% unnest() %>%
  mutate(StudyA = factor(StudyA,
                         levels = names(experiment_size)),
         StudyB = factor(StudyB,
                         levels = names(experiment_size)))


summary_concordance = jaccard_res_all %>% 
                      dplyr::filter(StudyA != StudyB) %>%
                      group_by(StudyA) %>%
                      summarise(mean(JaccardIx))

print("Summary concordance")
mean(summary_concordance[[2]])

##### Using top genes #####

study_genelist = lapply(METAheart, function(x){
  gene_list = dplyr::slice(x$HF_limma,1:500)[[1]]
  return(gene_list)
}) # get top 200 genes

# Make a matrix/data-frame with pairwaise comparisons
jaccard_res =  enframe(lapply(experiments, function(set_a){
  genes_a = study_genelist[[set_a]]
  j_ix_a = lapply(experiments, function(set_b){
    
    genes_b = study_genelist[[set_b]]
    
    #Jaccard Index
    j_ix = length(intersect(genes_a,genes_b))/length(union(genes_a,genes_b))
    
    return(j_ix)
  })
  
  j_ix_a = enframe(j_ix_a, "StudyB","JaccardIx") %>% 
    unnest()
  
  return(j_ix_a)
  
}), "StudyA") %>% unnest() %>%
  mutate(StudyA = factor(StudyA,
                         levels = names(experiment_size)),
         StudyB = factor(StudyB,
                         levels = names(experiment_size)))

# Data Manipulation for plots

jaccard_res = jaccard_res %>% spread(StudyA,JaccardIx) #First page excel

jaccard_res_mat = (as.matrix(jaccard_res[,-1]))
rownames(jaccard_res_mat) = jaccard_res[[1]]

jaccard_res_mat[upper.tri(jaccard_res_mat[names(experiment_size),
                                          names(experiment_size)])] = NA

jaccard_res_mat[jaccard_res_mat==1] = NA

jaccard_df = reshape2::melt(t(jaccard_res_mat),na.rm = T) %>% 
             mutate(Var2 = factor(as.character(Var2),
                                  levels = rev(names(experiment_size))))

saveRDS(jaccard_df, file = "data/figure_objects/jaccard_df.rds")

# Results reported

res_df = jaccard_df %>% dplyr::filter(Var1 != Var2)

print("Mean Jaccard Index")
mean(res_df$value)

#2. Transfer learning / disease score as classifier

t_matrix = get_all_limma(meta_list = METAheart,
                         limma_column = "t")

# Here we get AUC for pairwise classifiers
pairwise_200 = pairwise_ds(experiments = experiments,
                           meta_list = METAheart,
                           t_matrix = t_matrix,
                           ngenes = 500) #Second page excel

# Results from paper

pairwise_res = pairwise_200 %>%
  dplyr::filter(PredictorExperiment != 
                  PredictedExperiment)

print("Median AUROC")

median(pairwise_res$single)

print("Meann AUROC")

mean(pairwise_res$single)



saveRDS(pairwise_200[,c("PredictorExperiment", 
                        "PredictedExperiment", 
                        "single")], file = "data/figure_objects/pairwise_200.rds")

#3. Enrichment scores

study_deg_list_up = lapply(METAheart, function(x){
  deg = dplyr::slice(x$HF_limma,
                     1:500) %>%
    filter(t >0)
  return(deg[[1]])
})

study_deg_list_down = lapply(METAheart, function(x){
  deg = dplyr::slice(x$HF_limma,
                     1:500) %>%
    filter(t <0)
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

saveRDS(up_ES, 
        file = "data/figure_objects/up_ES.rds")


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

saveRDS(down_ES, 
        file = "data/figure_objects/down_ES.rds")

## Write supplementary tables for paper

WriteXLS(x = c("jaccard_res",
               "pairwise_200",
               "up_ES",
               "down_ES"), 
         ExcelFileName = "data/paper_sup/Reproducibility.xlsx",
         SheetNames = c("Jaccard_Index",
                        "AUC_disease_score",
                        "ES_upregulation",
                        "ES_downregulation")
)

## Correlations of reproducibility measures: 
## Here we want to know if high mean enrichment scores
## for a given study are correlated with their mean AUCs coming
## from the classification using the disease score

pairwise_200 = dplyr::filter(pairwise_200,
                             PredictorExperiment != PredictedExperiment) %>% 
               mutate(pairwise = paste(PredictorExperiment,
                                       PredictedExperiment,
                                       sep="_")) %>%
  dplyr::select(pairwise,single)

up_ES = dplyr::filter(up_ES,
                      DEG != Reference) %>%
              mutate(pairwise = paste(DEG,
                               Reference,
                               sep="_")) %>%
  dplyr::select(pairwise,ES)

median(up_ES$ES)
sum(up_ES$ES<0)


down_ES = dplyr::filter(down_ES,
                      DEG != Reference) %>%
         mutate(pairwise = paste(DEG,
                          Reference,
                          sep="_")) %>%
  dplyr::select(pairwise,ES)

median(down_ES$ES)
sum(down_ES$ES>0)

reprod_cor_data = left_join(pairwise_200,
          left_join(up_ES, down_ES, 
                    by = "pairwise"),
          by = "pairwise")

print("Correlation of ES with DS (up and down)")

cor.test(reprod_cor_data$single,reprod_cor_data$ES.x,method = "pearson")
cor.test(reprod_cor_data$single,reprod_cor_data$ES.y,method = "pearson")

#saveRDS(reprod_cor_data, 
#        file = "data/figure_objects/reprod_cor_data.rds")