# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of GSE1869

library(annotate)
library(hgu133a.db)
library(tidyverse)

#INCOMPLETE RAW DATA

#### Targets ####
GSE1869_targets = t(read.table("data_processing/raw/GSE1869/GSE1869_rawTargets.txt", sep ="\t",header = F,stringsAsFactors = F))
GSE1869_targets = as_tibble(GSE1869_targets[-1,c(2,12)])
colnames(GSE1869_targets) = c("Sample","Disease")
GSE1869_targets = mutate(GSE1869_targets, 
                         HeartFailure = ifelse(Disease=="Unused donor heart",
                                               "no","yes"))

GSE1869_targets = mutate(GSE1869_targets,
                          HTx = ifelse(grepl("Pre-LVAD",Disease),
                                       "no","yes"),
                          DCM = ifelse(grepl("nonischemic",Disease),
                                       "yes","no"))

GSE1869_targets[GSE1869_targets$HeartFailure == "no",
                "HTx"] = "yes"

#### counts ####

HF_expr_mat = read.table(file="data_processing/raw/GSE1869/GSE1869_rawcounts.txt",
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


## Incomplete RAW data set, complete processed data set ##
HF_expr_mat = HF_expr_mat[,GSE1869_targets$Sample]
GSE1869_counts = HF_expr_mat

save(GSE1869_counts, file = "data_processing/processed/GSE1869_counts.ro")
save(GSE1869_targets, file = "data_processing/processed/GSE1869_targets.ro")





