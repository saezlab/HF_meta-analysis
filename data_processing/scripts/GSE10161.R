# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE10161

library(limma)
library(oligo)
library(annotate)
library(hgu133a.db)
library(tidyverse)

#### Targets ####
GSE10161_targets = t(read.table("data_processing/raw/GSE10161/GSE10161_rawtargets.txt", 
                               sep ="\t",header = F,stringsAsFactors = F))
GSE10161_targets = as_tibble(GSE10161_targets[-1,c(2,8,11:20)])
colnames(GSE10161_targets)[c(1,2)] = c("Sample","Disease")
GSE10161_targets = mutate(GSE10161_targets, 
                         HeartFailure = ifelse(Disease=="control","no","yes"))

EF = c(GSE10161_targets$V8[grep("ef",GSE10161_targets$V8)],
       GSE10161_targets$V9[grep("ef",GSE10161_targets$V9)])

EF = as.numeric(gsub("ef: ","",EF))
age = as.numeric(gsub("age: ","",GSE10161_targets$V5))

GSE10161_targets = mutate(GSE10161_targets, EF = EF, age = age)

#### Count preproc ####
allfiles = list.files("data_processing/raw/GSE10161/GSE10161_RAW")
CELfiles = paste("data_processing/raw/GSE10161/GSE10161_RAW/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)
HF_exprnorm = rma(HF_expr)
HF_expr_mat = exprs(HF_exprnorm)

genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hgu133a.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),".C"),function(x) x[1]))

HF_expr_mat=HF_expr_mat[,GSE10161_targets$Sample]


##########################
##########################
GSE10161_counts = HF_expr_mat

save(GSE10161_counts, file = "data_processing/processed/GSE10161_counts.ro")
save(GSE10161_targets, file = "data_processing/processed/GSE10161_targets.ro")


