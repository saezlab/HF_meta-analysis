# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE42955

library(tidyverse)
library(limma)
library(oligo)
library(annotate)
library(hugene10sttranscriptcluster.db)

#### Targets ####
GSE42955_targets = t(read.table("data_processing/raw/GSE42955/GSE42955_rawtargets.txt", sep ="\t"
                                ,header = F,stringsAsFactors = F))
GSE42955_targets = as_tibble(GSE42955_targets[-1,c(2,8,10)])
colnames(GSE42955_targets) = c("Sample","Disease","Gender")
GSE42955_targets = mutate(GSE42955_targets, HeartFailure = ifelse(Disease=="Normal heart",
                                                                  "no","yes"))
GSE42955_targets = GSE42955_targets %>% mutate(Gender = unlist(lapply(strsplit(Gender,": "),
                                                                      function(x) x[2])))

GSE42955_targets = mutate(GSE42955_targets,
                          HTx = "yes",
                          DCM = ifelse(grepl("Dilated",Disease),
                                       "yes","no"))

#### Expression Processing ####
allfiles = list.files("data_processing/raw/GSE42955")
CELfiles = paste("data_processing/raw/GSE42955/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)

HF_exprnorm = rma(HF_expr,target="core")
HF_expr_mat = exprs(HF_exprnorm)

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hugene10sttranscriptcluster.db")

HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

#### Complete preproc ####
## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)

HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID


colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),"_"), 
                                      function(x) x[1]))

HF_expr_mat = HF_expr_mat[,GSE42955_targets$Sample]

GSE42955_counts = HF_expr_mat

save(GSE42955_counts, file = "data_processing/processed/GSE42955_counts.ro")
save(GSE42955_targets, file = "data_processing/processed/GSE42955_targets.ro")



