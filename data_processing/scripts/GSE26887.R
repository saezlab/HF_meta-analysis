# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE26887

library(tidyverse)
library(limma)
library(oligo)
library(annotate)
library(hugene10sttranscriptcluster.db)

#### Complete preproc ####
GSE26887_targets = t(read.table("data_processing/raw/GSE26887/GSE26887_rawTargets.txt", 
                                sep ="\t",header = F,stringsAsFactors = F))
GSE26887_targets = as_tibble(GSE26887_targets[-1,c(2,11,12,13)])
colnames(GSE26887_targets) = c("Sample","Gender","Age","Disease")

GSE26887_targets[2:ncol(GSE26887_targets)] = apply(GSE26887_targets[2:ncol(GSE26887_targets)],2,function(x){
                                              x = gsub(" ","",x)
                                              x= unlist(lapply((strsplit(x,":")), function(y) y[2]))
                                              return(x)
                                            })


GSE26887_targets = mutate(GSE26887_targets, HeartFailure = ifelse(Disease=="CONTROL" ,"no","yes"))

GSE26887_targets = GSE26887_targets %>% mutate(Gender = tolower(Gender),
                                               Age = as.numeric(gsub("years","",Age)))

GSE26887_targets = GSE26887_targets %>% 
                   mutate(HTx = ifelse(HeartFailure=="yes",
                                       "no","yes"),
                          DCM = ifelse(grepl("CONTROL",Disease),
                                       "no","yes"),
                          Diabetes = ifelse(grepl("NONDIABETIC",Disease) | 
                                                grepl("CONTROL",Disease),
                                              "no","yes"))



##############################################
allfiles = list.files("data_processing/raw/GSE26887")
CELfiles = paste("data_processing/raw/GSE26887/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
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
HF_expr_mat = HF_expr_mat[,GSE26887_targets$Sample]


#####################################
GSE26887_counts = HF_expr_mat

save(GSE26887_targets, file = "data_processing/processed/GSE26887_targets.ro")
save(GSE26887_counts, file = "data_processing/processed/GSE26887_counts.ro")
