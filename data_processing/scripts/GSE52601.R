# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE52601


library("illuminaHumanv4.db")
library(annotate)
library(tidyverse)

GSE52601_targets = t(read.table(file = "data_processing/raw/GSE52601/GSE52601_rawtargets.txt",
                               sep = "\t", ))[-1,c(2,8,11,12,13)]
colnames(GSE52601_targets) = c("Sample","Description","Disease","Age","Gender")

GSE52601_targets = as_tibble(data.frame(GSE52601_targets,stringsAsFactors = F))

GSE52601_targets = mutate(GSE52601_targets, 
       HeartFailure = ifelse(grepl("None",Disease),"no","yes"),
       Disease = unlist(lapply(strsplit(Disease,": "),
                               function(x) x[2])),
       Gender = unlist(lapply(strsplit(Gender,": "),
                              function(x) x[2])),
       Age = as.numeric(unlist(lapply(strsplit(Age,": "),
                                      function(x) x[2])))
)

# Filter fetal samples and after LVAD data
adult_GSE52601_targets = dplyr::filter(GSE52601_targets,
                                       (!grepl("fetal",Description))) %>%
  dplyr::filter(!(HeartFailure=="yes" & grepl("after",Description)))

# Filter HF samples for fetal targets
fetal_GSE52601_targets = dplyr::filter(GSE52601_targets,
                                      HeartFailure == "no") %>%
  mutate(fetal = ifelse(grepl("fetal",Description),
                        "yes","no"))

# Counts
GSE52601_matrix = read.table("data_processing/raw/GSE52601/GSE52601_counts.txt", sep = "\t", stringsAsFactors = F,
                            row.names = 1,header = T)

genesymbols = getSYMBOL(rownames(GSE52601_matrix), "illuminaHumanv4.db")

HF_expr_mat = as.matrix(GSE52601_matrix)
HF_expr_mat = HF_expr_mat[!is.na(genesymbols),] 
rownames(HF_expr_mat) = genesymbols[!is.na(genesymbols)]

## mean by same probeset
GENENAMES = rownames(HF_expr_mat)
HF_dataframe = data.frame(HF_expr_mat,stringsAsFactors = F) %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:ncol(HF_expr_mat)], 
                         by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

## Final check of order
adult_GSE52601_counts = HF_expr_mat[,adult_GSE52601_targets$Sample]
fetal_GSE52601_counts = HF_expr_mat[,fetal_GSE52601_targets$Sample]

# Here subset for what we need
save(adult_GSE52601_counts,file = "data_processing/processed/adult_GSE52601_counts.ro")
save(adult_GSE52601_targets,file = "data_processing/processed/adult_GSE52601_targets.ro")

save(fetal_GSE52601_counts,file = "data_processing/processed/fetal_GSE52601_counts.ro")
save(fetal_GSE52601_targets,file = "data_processing/processed/fetal_GSE52601_targets.ro")

