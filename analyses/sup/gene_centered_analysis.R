# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we test the associations of
#' PCs of the reduction of gene-centered matrices
#' with study labels

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(PLIER)
library(sjstats)

#0. Data processing

#Modifying meta analysis object to contain identical covariates and adding
#standardized data sets
METAheart = lapply(METAheart, function(x){
  x[["GEX_norm"]] = PLIER::rowNorm(x[["GEX"]]) #standarization of each data set (gene wise)
  return(x)
})

meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample, HeartFailure, Gender,
                Age, ExpID, HTx, DCM) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_"))

# Annotating RNAseq and Microarray
marrays = c("GSE76701","GSE57345","GSE42955",
            "GSE1869","GSE3585","GSE26887",
            "GSE5406","GSE16499")

load("./data/dictionaryIDs.ro")
new_ids = dictionary %>% dplyr::filter(GEO_ID %in% marrays)
marrays = new_ids$newID

meta_targets = meta_targets %>% mutate(Tech = ifelse(ExpID %in% marrays,
                                                     "microarray","rnaseq"))
# Generating unified data set (gene centered)
meta_gex_scale = get_complete_gex(meta_list = METAheart,
                                  complete_targets = meta_targets,
                                  gex_key = "GEX_norm")

#1. PCA in complete matrix
meta_gex_naomit_scale = meta_gex_scale[rowSums(is.na(meta_gex_scale))==0,]
pca_meta_scale = prcomp(t(meta_gex_naomit_scale[,meta_targets$grl_id]),
                        center = T,scale. = T)
pca_meta_sum_scale = summary(pca_meta_scale)
pca_plot_df_scale = meta_targets %>% mutate("PC1" = pca_meta_scale$x[,1],
                                            "PC2" = pca_meta_scale$x[,2])

pca_meta_sum_scale[["plot_df"]] = pca_plot_df_scale

saveRDS(pca_meta_sum_scale, 
        file = "data/figure_objects/gcentered_PCA_sum.rds")

#2. Association of PCs with HF
pca_anova_res = apply((t(pca_meta_scale$x)), 1, function(x, targets){
  pc_i = x
  #factor_a_vect = factor(targets[["ExpID"]])
  factor_a_vect = factor(targets[["HeartFailure"]])
  gene_aov = aov(pc_i ~ factor_a_vect)
  aov_stats = anova_stats(gene_aov) 
  
},targets = meta_targets)  %>% bind_rows(.id = "PC") %>% as_tibble() %>% 
  group_by(PC) %>% gather(stats,value,-(PC:term))  %>% spread(term,value) %>%
  ungroup()

pcs_hf = pca_anova_res %>% filter(stats == "p.value" & factor_a_vect < 0.05)

pcs_hf = pca_meta_sum_scale$importance[,pcs_hf$PC]

total_prop_hf = sum(pcs_hf[2,])

print(total_prop_hf)

#3. Association of PCs with Study
pca_anova_res = apply((t(pca_meta_scale$x)), 1, function(x, targets){
  pc_i = x
  factor_a_vect = factor(targets[["ExpID"]])
  gene_aov = aov(pc_i ~ factor_a_vect)
  aov_stats = anova_stats(gene_aov) 
  
},targets = meta_targets)  %>% bind_rows(.id = "PC") %>% as_tibble() %>% 
  group_by(PC) %>% gather(stats,value,-(PC:term))  %>% spread(term,value) %>%
  ungroup()

pcs_study = pca_anova_res %>% filter(stats == "p.value" & factor_a_vect < 0.05)

pcs_study = pca_meta_sum_scale$importance[,pcs_study$PC]

total_prop_study = sum(pcs_study[2,])

print(total_prop_study)

#4. Creation of table for plot
summary_table = cbind(c("Heart Failure","study"),rbind(c(ncol(pcs_hf), total_prop_hf),
                      c(ncol(pcs_study), total_prop_study)))

colnames(summary_table) = c("factor","Total PCs","Cumulative Proportion \n of Variance")
rownames(summary_table) = c("Heart Failure", "Study")

saveRDS(summary_table, 
        file = "data/figure_objects/gcentered_PCs_sum.rds")



































