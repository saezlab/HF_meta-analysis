# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we make the list object that includes
#' all especial datasets not included in the meta-analysis
#' but used in the validation analysis

library(tidyr)
library(dplyr)

#### VALIDATION STUDIES ####
ValidationStudies = list()

##GSE84796
load(file = "data_processing/processed/GSE84796_counts.ro")
load(file = "data_processing/processed/GSE84796_targets.ro")

ValidationStudies[["GSE84796"]] = list("GEX"= GSE84796_matrix,
                                       "TARGETS" = GSE84796_targets)

##GSE4172
load("data_processing/processed/GSE4172_counts.ro")
load("data_processing/processed/GSE4172_targets.ro")

ValidationStudies[["GSE4172"]] = list("GEX"= GSE4172_counts,
                                      "TARGETS" = GSE4172_targets)

#GSE9800
load("data_processing/processed/GSE9800_counts.ro")
load("data_processing/processed/GSE9800_targets.ro")

ValidationStudies[["GSE9800"]] = list("GEX"= GSE9800_counts,
                                      "TARGETS" = GSE9800_targets)

#GSE10161
load("data_processing/processed/GSE10161_counts.ro")
load("data_processing/processed/GSE10161_targets.ro")
ValidationStudies[["GSE10161"]] = list("GEX"= GSE10161_counts,
                                       "TARGETS" = GSE10161_targets)

#GSE52601
load("data_processing/processed/adult_GSE52601_counts.ro")
load("data_processing/processed/adult_GSE52601_targets.ro")

ValidationStudies[["GSE52601"]] = list("GEX"= adult_GSE52601_counts,
                                       "TARGETS" = adult_GSE52601_targets)

#GSE3586 - from processed - Gender/Age *
load("data_processing/processed/GSE3586_matrix.ro")
load("data_processing/processed/GSE3586_targets.ro") 

GSE3586_target = GSE3586_targets %>%
  dplyr::filter(Sample %in% colnames(GSE3586_matrix))

ValidationStudies[["GSE3586"]] = list("GEX" = GSE3586_matrix,
                                      "TARGETS" = GSE3586_target)

# Creating object

ValidationStudies = lapply(ValidationStudies,function(x){
  x$GEX = x$GEX[,x$TARGETS$Sample]
  return(x)
})

saveRDS(ValidationStudies, file = "data/external_METAheart.rds")