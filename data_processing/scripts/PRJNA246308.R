# Author: Jan Lanzer, 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing Pipeline for Study PRJNA246308
# License: MIT 

library(Biobase)
library(tidyverse)
library(limma)
library(edgeR)

# 1)  Read sample information, process, and save as target file. 
PRJNA246308_target = read.table("HGEX_data/Study_Processing/Metainformation/PRJNA246308.txt", 
                                header = TRUE, 
                                stringsAsFactors = F)

# add information manually curated from PRJNA246308_targetfile
Gender = c("female","male", "male", "male","female")
Age = c(57,NA, 52, 68, 59)

# edit information and format columns
PRJNA246308_target = as_tibble(cbind(PRJNA246308_target,Gender,Age)) %>% 
  mutate(Gender = as.character(Gender)) %>% 
  dplyr::rename(Sample = sample_alias) %>% 
  dplyr::rename(Disease = sample_title) %>% 
  mutate(HeartFailure = "yes")

PRJNA246308_target = PRJNA246308_target%>% select(Sample, Disease, HeartFailure, Gender, Age) %>%
  mutate(DCM = "yes") %>% mutate(HTx = "yes")

PRJNA246308_target$HeartFailure[str_detect(PRJNA246308_target$Disease, "NF")] = "no"
PRJNA246308_target$DCM[str_detect(PRJNA246308_target$Disease, "NF")] = "no"

#save r-object
  save(PRJNA246308_target, file = "HGEX_data/Study_Processing/PRJNA246308_target.ro")


# 2) Read gene expression data, filter & normalize, and save as count file. 
# GEX table
df = read.table("HGEX_data/Study_Processing/GEXraw/PRJNA246308-rawdata.txt")
colnames(df)= PRJNA246308_target$Sample

#### Filter & Normalize
#create DGE class object
group <- PRJNA246308_target$HeartFailure
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

### save R-objects
PRJNA246308_count= v$E
save(PRJNA246308_count, file = "HGEX_data/Study_Processing/PRJNA246308_count.ro")



