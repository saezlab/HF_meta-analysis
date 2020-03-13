# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE16499

library(tidyverse)
library(limma)
library(oligo)
library(annotate)
library(huex10sttranscriptcluster.db)

#### Target preproc ####
GSE16499_targets = t(read.table("data_processing/raw/GSE16499/GSE16499_rawtargets.txt", 
                                sep ="\t",header = F,stringsAsFactors = F))
GSE16499_targets = as_tibble(GSE16499_targets[-1,c(2,1)])
colnames(GSE16499_targets) = c("Sample","Disease")
GSE16499_targets$DCM = "no"
GSE16499_targets$HTx = "yes"
GSE16499_targets$HeartFailure = ifelse(grepl("Normal",GSE16499_targets$Disease),
                                       "no","yes")

#### Raw files ####
allfiles = list.files("data_processing/raw/GSE16499")
CELfiles = paste("data_processing/raw/GSE16499/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)

HF_exprnorm = rma(HF_expr,target="core")
HF_expr_mat = exprs(HF_exprnorm)
genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "huex10sttranscriptcluster.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)

HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),".CEL."), 
                                      function(x) x[1]))

HF_expr_mat = HF_expr_mat[,GSE16499_targets$Sample]

##
GSE16499_counts = HF_expr_mat

save(GSE16499_targets, file = "data_processing/processed/GSE16499_targets.ro")
save(GSE16499_counts, file = "data_processing/processed/GSE16499_counts.ro")




