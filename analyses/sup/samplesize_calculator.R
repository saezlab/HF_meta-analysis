# MIT License

# Copyright (c) [2020] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: claculate sample sizes and gene coverage

library(tidyverse)

META  = readRDS("HGEX_data/METAheart.rds")

ct= c()
dcm= c()
icm = c()
genes = c()

for (study in names(META)){
  ct= c(ct,dim(META[[study]]$TARGETS %>% filter(HeartFailure == "no"))[1])
  dcm = c(dcm,dim(META[[study]]$TARGETS %>% filter(DCM == "yes"))[1])
  icm = c(icm,dim(META[[study]]$TARGETS %>% filter(DCM == "no") %>% filter(HeartFailure == "yes"))[1])
  genes = c(genes, dim(META[[study]]$GEX)[1])
  }

sample.size = data.frame("study"= as.character(names(META)), "CT"= ct, "DCM"= dcm, "ICM"= icm, "genes" = genes) %>% 
  mutate(total = rowSums(.[,2:4])) %>%
  arrange(desc(total))

saveRDS(sample.size, file ="HGEX_data/clinical_description/sample_sizes.rds")




