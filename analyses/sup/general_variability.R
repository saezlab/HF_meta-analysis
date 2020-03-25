# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform PCAs of merged data,
#' z-transformed HF samples and generate associations
#' of pcs with known covariates

#Inputs and dependencies
source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(PLIER)
library(Rtsne)
library(sjstats)

# All samples' information
meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample, HeartFailure, Gender,
                Age, ExpID, HTx, DCM) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_"))

# A single gene expression matrix with all samples
meta_gex = get_complete_gex(meta_list = METAheart,
                            complete_targets = meta_targets,
                            gex_key = "GEX")

# Annotating RNAseq and Microarray

marrays = c("GSE76701","GSE57345","GSE42955",
            "GSE1869","GSE3585","GSE26887",
            "GSE5406","GSE16499")

load("data/dictionaryIDs.ro")
new_ids = dictionary %>% dplyr::filter(GEO_ID %in% marrays)
marrays = new_ids$newID

meta_targets = meta_targets %>% mutate(Tech = ifelse(ExpID %in% marrays,
                                                     "microarray","rnaseq"))

# 1. PCA of all data (complete genes)

meta_gex_naomit = meta_gex[rowSums(is.na(meta_gex))==0,]
dim(meta_gex_naomit)
pca_meta = prcomp(t(meta_gex_naomit[,meta_targets$grl_id]),center = T,scale. = T)
pca_meta_sum = summary(pca_meta) #This is the object to plot

pca_meta_sum[["plot_df"]] = meta_targets %>% mutate("PC1" = pca_meta$x[,1],
                                                     "PC2" = pca_meta$x[,2])

saveRDS(pca_meta_sum, file = "data/figure_objects/pca_meta_summary.rds")

# 2. Matrix of z-scores + PCA

METAheart = lapply(METAheart, function(x){
  
  targets = x[["TARGETS"]]
  gex = x[["GEX"]]
  
  hf_gex = gex[,(targets %>% 
                 dplyr::filter(HeartFailure == "yes"))$Sample]

  healthy_gex = gex[,(targets %>% 
                   dplyr::filter(HeartFailure == "no"))$Sample]
  
  ref_mean = rowMeans(healthy_gex)
  
  ref_sd = apply(healthy_gex,1,sd)
  
  x[["Zmat"]] = (hf_gex - ref_mean)/ref_sd
  
  return(x)
  
})

z_targets = meta_targets %>% filter(HeartFailure == "yes")

meta_gex_z = get_complete_gex(meta_list = METAheart,
                            complete_targets = z_targets,
                            gex_key = "Zmat")

meta_gex_z_naomit = meta_gex_z[rowSums(is.na(meta_gex_z))==0,]

pca_meta_z = prcomp(t(meta_gex_z_naomit[,z_targets$grl_id]),
                    center = T,scale. = T)

pca_meta_z_sum = summary(pca_meta_z) #This is the object to plot

pca_meta_z_sum[["plot_df"]] = z_targets %>% mutate("PC1" = pca_meta_z$x[,1],
                                                     "PC2" = pca_meta_z$x[,2])

saveRDS(pca_meta_z_sum, file = "data/figure_objects/pca_meta_summary_z.rds")

# 3. Fitting a linear model to each PC of z-transformed hf samples
# to find associations with study

pcs_study = run_anovastats_single(numeric_matrix = t(pca_meta_z_sum$x),
                                  targets = z_targets,
                                  factor_a = "ExpID",
                                  pval = 0.05)

# Map to proportion explained

pcs_study = mutate(pcs_study,
                   prop_var = pca_meta_z_sum$importance[2,pcs_study$PC]) %>%
  dplyr::arrange(desc(prop_var)) %>% 
  dplyr::select(PC, factor_a_vect, prop_var) %>% 
  dplyr::mutate(factor_a_vect = as.character(factor_a_vect))

total_row = c("TOTAL","",sum(pcs_study$prop_var))
pcs_study = pcs_study %>% mutate(prop_var = as.character(prop_var))
names(total_row) = colnames(pcs_study) = c("PCs","p_value","Prop. Var")

pcs_study = bind_rows(pcs_study, total_row)

print(pcs_study,n=100)

# 4. Fitting a linear model to each PC of z-transformed hf samples
# to find associations with DCM vs nDCM

pcs_dcm = run_anovastats_single(numeric_matrix = t(pca_meta_z_sum$x),
                                targets = z_targets,
                                factor_a = "DCM",
                                pval = 0.05)

# Map to proportion explained

pcs_dcm = mutate(pcs_dcm,
                   prop_var = pca_meta_z_sum$importance[2,pcs_dcm$PC]) %>%
  dplyr::arrange(desc(prop_var)) %>% 
  dplyr::select(PC, factor_a_vect, prop_var) %>% 
  dplyr::mutate(factor_a_vect = as.character(factor_a_vect))

total_row = c("TOTAL","",sum(pcs_dcm$prop_var))
pcs_dcm = pcs_dcm %>% mutate(prop_var = as.character(prop_var))
names(total_row) = colnames(pcs_dcm) = c("PCs","p_value","Prop. Var")

pcs_dcm = bind_rows(pcs_dcm, total_row)

# 5. Fitting a linear model to each PC of z-transformed hf samples
# to find associations with HTx

pcs_htx = run_anovastats_single(numeric_matrix = t(pca_meta_z_sum$x),
                                targets = z_targets,
                                factor_a = "HTx",
                                pval = 0.05)

# Map to proportion explained

pcs_htx = mutate(pcs_htx,
                 prop_var = pca_meta_z_sum$importance[2,pcs_htx$PC]) %>%
  dplyr::arrange(desc(prop_var)) %>% 
  dplyr::select(PC, factor_a_vect, prop_var) %>% 
  dplyr::mutate(factor_a_vect = as.character(factor_a_vect))

total_row = c("TOTAL","",sum(pcs_htx$prop_var))
pcs_htx = pcs_htx %>% mutate(prop_var = as.character(prop_var))
names(total_row) = colnames(pcs_htx) = c("PCs","p_value","Prop. Var")

pcs_htx = bind_rows(pcs_htx, total_row)

# 5. t-SNE of z-transformed data

## Rtsne function may take some minutes to complete...
set.seed(9)  
tsne_model = Rtsne(t(meta_gex_z_naomit[,z_targets$grl_id]),
                     pca=TRUE, dims=2)

## getting the two dimension matrix
d_tsne = as.data.frame(tsne_model$Y) 

tsne_plotdf = z_targets %>% mutate("tSNE1" = d_tsne[,1],
                                   "tSNE2" = d_tsne[,2])

saveRDS(tsne_plotdf, file = "data/figure_objects/tsne_z.rds")





















