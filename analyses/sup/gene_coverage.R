# MIT License

# Copyright (c) [2020] [Ricardo O. Ramirez, Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: Comparing gene coverage of studies included in meta analysis 
# by calculating Jaccard Index of all genes


library(tidyverse)

#load data
METAheart= readRDS("HGEX_data/METAheart.rds")
experiments = names(METAheart)
names(experiments) = experiments

# Jaccard index of all genes
# 1. Jaccaard Index
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

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

# mean Jaccardindex, added to figure manually
print("Summary concordance") 
mean(summary_concordance[[2]])

## continue data transformation
jaccard_res_all = jaccard_res_all %>% 
  mutate(JaccardIx= ifelse(StudyA==StudyB, NA,JaccardIx))%>% #add NAs to same studies
  spread(StudyA,JaccardIx) 

jaccard_res_mat = (as.matrix(jaccard_res_all[,-1]))
rownames(jaccard_res_mat) = jaccard_res_all[[1]]

jaccard_res_mat[upper.tri(jaccard_res_mat[names(experiment_size),
                                          names(experiment_size)])] = NA

jaccard_df = reshape2::melt(t(jaccard_res_mat),na.rm = T) %>% 
  mutate(Var2 = factor(as.character(Var2),
                       levels = rev(names(experiment_size))))

#save file for plotting
saveRDS(jaccard_df, file = "HGEX_data/figure_objects/jaccard_df.rds")


