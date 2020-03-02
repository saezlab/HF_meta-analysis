# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE57345

library(tidyverse)
library(limma)
library(annotate)
library(hugene11sttranscriptcluster.db)
library(oligo)

# GSE57345
rawtargets = read.table("data_processing/raw/GSE57345/TargetAnnotation.txt",sep = "\t", stringsAsFactors = F)
rawtargets = t(rawtargets)
colnames(rawtargets) = rawtargets[1,]
rawtargets = rawtargets[-1,]

Targets = rawtargets[,c(2,11,12,13,14)]
Targets = as_tibble(Targets)

Targets = mutate(Targets, HeartFailure = unlist(lapply(strsplit(Targets[[2]],": "),
                                       function(x) x[2])),
                 Disease = unlist(lapply(strsplit(Targets[[3]],": "),
                                  function(x) x[2])),
                 Gender = unlist(lapply(strsplit(Targets[[4]],": "),
                                 function(x) x[2])),
                 Age = as.numeric(unlist(lapply(strsplit(Targets[[5]],": "),
                                     function(x) x[2])))
                 )

Targets = mutate(Targets, Disease = gsub("-","_",gsub(" ","_", Targets$Disease)))
colnames(Targets)[1] = "Sample"

GSE57345_targets = Targets[,c("Sample", "HeartFailure",
                              "Disease","Gender","Age")]

GSE57345_targets = mutate(GSE57345_targets,
                          HTx = "yes",
                          DCM = ifelse(grepl("dilated",Disease),
                                        "yes","no"))

save(GSE57345_targets, file = "data_processing/processed/GSE57345_targets.ro")

#####################################################################3

allfiles = list.files("data_processing/raw/GSE57345")
CELfiles = paste("data_processing/raw/GSE57345/",sort(allfiles[grep(".CEL",allfiles)]),sep="")
HF_expr = read.celfiles(CELfiles)
HF_exprnorm = rma(HF_expr,target="core")
HF_expr_mat = exprs(HF_exprnorm)

#Names
genesymbols = getSYMBOL(as.character(rownames(HF_expr_mat)), "hugene11sttranscriptcluster.db")
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),]
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

#### Complete preproc ####
## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], 
                         by = list(ID = HF_dataframe$ID), 
                         FUN = "mean", na.rm = T)

HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID


colnames(HF_expr_mat) = unlist(lapply(strsplit(colnames(HF_expr_mat),"_"), 
                                      function(x) x[1]))

HF_expr_mat = HF_expr_mat[,GSE57345_targets$Sample]

HF_mat = HF_expr_mat

GSE57345_counts = HF_expr_mat
        
save(GSE57345_counts, file = "data_processing/processed/GSE57345_counts.ro")
