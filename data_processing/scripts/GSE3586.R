# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE3586

library(biomaRt)
library(tidyverse)
library(limma)

GSE3586_targets = t(read.table("data_processing/raw/GSE3586/GSE3586_targets.txt", sep="\t",stringsAsFactors = F))
GSE3586_targets = data.frame(GSE3586_targets[-1,c(2,10)],stringsAsFactors = F)
colnames(GSE3586_targets) = c("Sample","Disease")
GSE3586_targets = mutate(GSE3586_targets,
                         HeartFailure = ifelse(Disease == "Non failing", 
                                               "no","yes")) %>% as_tibble()

####################################################################################################3

GSE3586_matrix = read.table("data_processing/raw/GSE3586/GSE3586_matrix.txt", sep = "\t", stringsAsFactors = F,
                            row.names = 1,header = T)
array_ids = read.table("data_processing/raw/GSE3586/array_annotation.txt", 
                       sep="\t",stringsAsFactors = F,header = T,comment.char = "#")[,c(1,14)]

array_ids = filter(array_ids, ID %in% rownames(GSE3586_matrix) & GENE_SYMBOL!="")
array_ids = filter(array_ids,GENE_SYMBOL!="")
GSE3586_matrix = GSE3586_matrix[array_ids$ID,]

HF_expr_mat = GSE3586_matrix

## mean by same probeset
GENENAMES = array_ids$GENE_SYMBOL
HF_dataframe = mutate(HF_expr_mat, ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], 
                         by = list(ID = HF_dataframe$ID), 
                         FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID
HF_expr_mat = HF_expr_mat[,GSE3586_targets$Sample]

####################################################################################################3

GSE3586_matrix = HF_expr_mat

save(GSE3586_matrix,file = "data_processing/processed/GSE3586_matrix.ro")
save(GSE3586_targets,file = "data_processing/processed/GSE3586_targets.ro")
















