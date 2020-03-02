# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE9800

library(tidyverse)
library(limma)

# Targets
GSE9800_targets = t(read.table("data_processing/raw/GSE9800/GSE9800_rawtargets.txt",
                               sep ="\t",header = F,stringsAsFactors = F))
GSE9800_targets = as_tibble(GSE9800_targets[-1,])
colnames(GSE9800_targets) = c("Disease","Sample")
GSE9800_targets = mutate(GSE9800_targets,Disease = gsub(" \\(.+\\)","",Disease),
                         HeartFailure = ifelse(grepl("patients",Disease), "yes","no"))

# Counts
GSE9800_counts = as.matrix((read.table("data_processing/raw/GSE9800/GSE9800_rawcounts.txt",
                                       sep ="\t",header = T,stringsAsFactors = F,
                             row.names = 1)))

GSE9800_annotations = (read.table("data_processing/raw/GSE9800/ANNOTATION.txt",
                                  sep ="\t",header = T,stringsAsFactors = F) %>% filter(GENE_SYMBOL!= ""))
GSE9800_counts =  GSE9800_counts[as.character(GSE9800_annotations$ID),]
rownames(GSE9800_counts) = GSE9800_annotations$GENE_SYMBOL


GSE9800_counts = GSE9800_counts[,-ncol(GSE9800_counts)]

######################################################

GSE9800_targets = filter(GSE9800_targets, Sample %in% colnames(GSE9800_counts))
GSE9800_counts = GSE9800_counts[,GSE9800_targets$Sample]

######################################################
HF_expr_mat = GSE9800_counts
#### Complete preproc ####
## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], 
                         by = list(ID = HF_dataframe$ID), 
                         FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

GSE9800_counts = HF_expr_mat
GSE9800_counts = GSE9800_counts[,GSE9800_targets$Sample]

#######################################################
save(GSE9800_counts,file = "data_processing/processed/GSE9800_counts.ro")
save(GSE9800_targets,file = "data_processing/processed/GSE9800_targets.ro")




















