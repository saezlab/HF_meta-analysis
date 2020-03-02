# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE4172

library(annotate)
library(hgu133a2.db)
library(tidyverse)
library(limma)

GSE4172_targets = t(read.table("data_processing/raw/GSE4172/GSE4172_targets.txt", 
                               sep ="\t",header = F,stringsAsFactors = F))
GSE4172_targets = as_tibble(GSE4172_targets[-1,c(1,7,9,10,11,12,13)])
colnames(GSE4172_targets) = c("Sample","Disease","Age","Gender","EjectionFraction","Diameter","Inflammation")

GSE4172_targets[,3:ncol(GSE4172_targets)] = apply(GSE4172_targets[,3:ncol(GSE4172_targets)],2, function(x){
  
  x = gsub(" ","",x)
  x = unlist(lapply(strsplit(x,":"), function(y) y[2]))
  
})


GSE4172_targets = mutate(GSE4172_targets, HeartFailure = ifelse(Disease == "healthy control", "no", "yes"))

###########
HF_expr_mat = read.table(file="data_processing/raw/GSE4172/GSE4172_counts.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
rownames(HF_expr_mat) = HF_expr_mat[,1]
HF_expr_mat = as.matrix(HF_expr_mat[,-1])

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hgu133a2.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset

GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)

HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID
HF_expr_mat = log2(HF_expr_mat)

HF_expr_mat = HF_expr_mat[,GSE4172_targets$Sample]

###########
GSE4172_counts = HF_expr_mat

save(GSE4172_counts, file = "data_processing/processed/GSE4172_counts.ro")
save(GSE4172_targets, file = "data_processing/processed/GSE4172_targets.ro")








