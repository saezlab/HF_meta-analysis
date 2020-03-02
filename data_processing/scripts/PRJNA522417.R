# Author: Jan Lanzer, 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing Pipeline for Study PRJNA522417
# License: MIT 

library(Biobase)
library(tidyverse)
library(limma)
library(edgeR)

# 1)  Read sample information, process, and save as target file. 
meta =as_tibble(read.table("HGEX_data/Study_Processing/Metainformation/GSE126569_information_edit.txt", header = TRUE, sep = '\t'))
meta2 = as_tibble(read.table("HGEX_data/Study_Processing/Metainformation/PRJNA522417.txt", header =TRUE, sep = '\t', stringsAsFactors = F)) %>% 
  select(experiment_accession,run_accession) %>% 
  dplyr::rename(SRA =experiment_accession)

# Update meta with HeartFailure and Developmental Information
meta = meta %>% mutate(HeartFailure = c(rep("no", 18), rep("yes", 15), rep("no", 5))) %>% 
  mutate(fetal= c(rep("no", 33), rep("yes",5))) %>%
  select(HeartFailure, fetal, SRA) %>%
  mutate(SRA = as.character(SRA)) %>%
  inner_join(meta2, by= "SRA") %>%
  rename(Sample = run_accession)

# Target file subsetting
PRJNA522417_target_fetal = meta %>% filter(meta$fetal== "yes")
PRJNA522417_target_adult=  meta %>% filter(meta$fetal== "no")

PRJNA522417_target_adult = PRJNA522417_target_adult %>% select(-SRA) %>%
  mutate(DCM = HeartFailure) %>%
  mutate(HTx = "yes") 

# Save target files
save(PRJNA522417_target_adult, file = "HGEX_data/Study_Processing/PRJNA522417_target_adult.ro")
save(PRJNA522417_target_fetal, file = "HGEX_data/Study_Processing/PRJNA522417_target_fetal.ro")


# 2) Read gene expression data, filter & normalize, and save as count file. 
# read GEX table
df = read.table("HGEX_data/Study_Processing/GEXraw/PRJNA522417-rawdata.txt")

#create DGE class object
group <- c(rep(1,18),rep(2,15),rep(3,5))
dge<- DGEList(counts=df, group=group)

# filter low expressed gene. Filtering keeps genes that have count-per-million (CPM) above k in n
# samples, where k is determined by min.count and by the sample library sizes and n is determined by the design matrix.
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization TMM. TMM is the recommended for most RNA-Seq data where the majority (more than half) 
#of the genes are believed not differentially expressed between any pair of the samples.
dge <- calcNormFactors(dge)

# use limma voom
v <- voom(dge, plot=TRUE)

# save gene count table
PRJNA522417_count= v$E

#count file subsetting
PRJNA522417_count_fetal= subset(PRJNA522417_count, select = c(PRJNA522417_target_fetal$Sample))
PRJNA522417_count_adult= PRJNA522417_count[,PRJNA522417_target_adult$Sample]

# save R-object
save(PRJNA522417_count_adult, file = "HGEX_data/Study_Processing/PRJNA522417_count_adult.ro")
save(PRJNA522417_count_fetal, file = "HGEX_data/Study_Processing/PRJNA522417_count_fetal.ro")
