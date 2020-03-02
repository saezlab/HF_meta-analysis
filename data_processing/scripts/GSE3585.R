# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE3585

library(limma)
library(oligo)
library(annotate)
library(hgu133a.db)
library(tidyverse)

#### Targets ####
GSE3585_targets = t(read.table("data_processing/raw/GSE3585/GSE3585_rawTargets.txt", 
                               sep ="\t",header = F,stringsAsFactors = F))
GSE3585_targets = as_tibble(GSE3585_targets[-1,c(2,10)])
colnames(GSE3585_targets) = c("Sample","Disease")
GSE3585_targets = mutate(GSE3585_targets, 
                         HeartFailure = ifelse(Disease=="Non failing","no","yes"))

#### Count preproc ####
allfiles = list.files("data_processing/raw/GSE3585")
CELfiles = paste("data_processing/raw//GSE3585/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)
HF_exprnorm = rma(HF_expr)
HF_expr_mat = exprs(HF_exprnorm)

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hgu133a.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], 
                         by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),".C"),function(x) x[1]))

HF_expr_mat=HF_expr_mat[,GSE3585_targets$Sample]

##########################
GSE3585_counts = HF_expr_mat
save(GSE3585_counts, file = "data_processing/processed/GSE3585_counts.ro")
#save(GSE3585_targets, file = "HGEX_data/src/GSE3585/preproc/GSE3585_targets.ro") #We better used the paper supplemental table




