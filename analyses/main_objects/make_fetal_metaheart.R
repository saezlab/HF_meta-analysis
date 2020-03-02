# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we make the list object that includes
#' all fetal datasets not included in the meta-analysis
#' but used in the validation analysis

library(tidyverse)
library(limma)

fetal_METAheart = list()

##############

load(file = "data_processing/processed/fetal_GSE52601_counts.ro")
load(file = "data_processing/processed/fetal_GSE52601_targets.ro")

fetal_METAheart[["GSE52601"]] = list("GEX" = fetal_GSE52601_counts,
                                     "TARGETS" = fetal_GSE52601_targets)

############
load(file = "data_processing/processed/PRJNA522417_target_fetal.ro")
load(file = "data_processing/processed/PRJNA522417_count_fetal.ro")
load(file = "data_processing/processed/PRJNA522417_count_adult.ro")
load(file = "data_processing/processed/PRJNA522417_target_adult.ro")

adult_healthy = filter(PRJNA522417_target_adult,HeartFailure == "no") %>% 
  dplyr::select(Sample, HeartFailure,fetal)

adult_count = PRJNA522417_count_adult[,adult_healthy$Sample]

fetal_targets = PRJNA522417_target_fetal %>% dplyr::select(Sample, HeartFailure,fetal)

PRJNA522417_fetal_targets = bind_rows(adult_healthy,fetal_targets)

PRJNA522417_fetal_counts = cbind(PRJNA522417_count_adult[,adult_healthy$Sample],
                                 PRJNA522417_count_fetal[,fetal_targets$Sample])

fetal_METAheart[["PRJNA522417"]] = list("GEX" = PRJNA522417_fetal_counts,
                                        "TARGETS" = PRJNA522417_fetal_targets)

saveRDS(fetal_METAheart, file =  "data/fetal_METAheart.rds")