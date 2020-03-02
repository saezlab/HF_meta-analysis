# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE5406

library(oligo)
library(annotate)
library(hgu133a.db)
library(tidyverse)

#### Targets ####
GSE5406_targets = t(read.table("data_processing/raw/GSE5406/GSE5406_rawtargets.txt", sep ="\t",header = F,stringsAsFactors = F))
colnames(GSE5406_targets) = c("Disease","Sample")
GSE5406_targets = as_tibble(GSE5406_targets[-1,c(2,1)])
GSE5406_targets = mutate(GSE5406_targets, HeartFailure = ifelse(grepl("NonFailing",Disease),"no","yes"))

GSE5406_targets = mutate(GSE5406_targets,
                          HTx = "yes",
                          DCM = ifelse(grepl("Idiopathic",Disease),
                                       "yes","no"))

#### ExprMat ####
HF_expr_mat = read.table(file="data_processing/raw/GSE5406/GSE5406_rawcounts.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
rownames(HF_expr_mat) = HF_expr_mat[,1]
HF_expr_mat = as.matrix(HF_expr_mat[,-1])

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hgu133a.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)

HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

HF_expr_mat = HF_expr_mat[,GSE5406_targets$Sample]

GSE5406_counts = HF_expr_mat

save(GSE5406_counts, file = "data_processing/processed/GSE5406_counts.ro")
save(GSE5406_targets, file = "data_processing/processed/GSE5406_targets.ro")


