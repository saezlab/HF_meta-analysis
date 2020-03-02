# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE76701

library(limma)
library(oligo)
library(annotate)
library(hgu133plus2.db)
library(tidyverse)

#### Targets ####
GSE76701_targets = t(read.table("data_processing/raw/GSE76701/GSE76701_rawtargets.txt", 
                               sep ="\t",header = F,stringsAsFactors = F))
GSE76701_targets = as_tibble(GSE76701_targets[-1,c(1,7)])
colnames(GSE76701_targets) = c("Sample","Disease")
GSE76701_targets = mutate(GSE76701_targets, 
                         HeartFailure = ifelse(grepl("Failing",Disease),"yes","no"),
                         DCM = "no")

#### Count preproc ####
allfiles = list.files("data_processing/raw/GSE76701")
CELfiles = paste("data_processing/raw/GSE76701/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)
HF_exprnorm = rma(HF_expr)
HF_expr_mat = exprs(HF_exprnorm)

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hgu133plus2.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),"_"),function(x) x[1]))

HF_expr_mat=HF_expr_mat[,GSE76701_targets$Sample]
plotMDS(HF_expr_mat,labels = GSE76701_targets$HeartFailure)

GSE76701_counts = HF_expr_mat

save(GSE76701_targets, file = "data_processing/processed/GSE76701_targets.ro")
save(GSE76701_counts, file = "data_processing/processed/GSE76701_counts.ro")

