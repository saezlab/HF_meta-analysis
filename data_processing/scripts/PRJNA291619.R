# Author: Jan Lanzer, 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing Pipeline for Study PRJNA291619
# License: MIT 

library(Biobase)
library(dplyr)
library(limma)
library(edgeR)
library(tidyverse)

# 1)  Read sample information, process, and save as target file. 
PRJNA291619_target = read.table("HGEX_data/Study_Processing/Metainformation/PRJNA291619.txt", 
                                header = TRUE,
                                stringsAsFactors = F)

# format columns
PRJNA291619_target = PRJNA291619_target %>%
  dplyr::rename(Sample = sample_alias) %>% 
  mutate(HeartFailure = "yes") %>% 
  dplyr::rename(Disease = sample_title) %>% 
  select(-study_accession,-run_accession) %>% 
  mutate(DCM = "no") %>% 
  mutate(HTx ="yes")

PRJNA291619_target$HeartFailure[str_detect(PRJNA291619_target$Disease, "Control")] = "no"
PRJNA291619_target$DCM[str_detect(PRJNA291619_target$Disease, "dilated")] = "yes"

# remove RCM samples
PRJNA291619_target = PRJNA291619_target[-c(5,6),]

# save R object
save(PRJNA291619_target, file = "HGEX_data/Study_Processing/PRJNA291619_target.ro")


# 2) Read gene expression data, filter & normalize, and save as count file. 
# read GEX table
df = read.table("HGEX_data/Study_Processing/GEXraw/PRJNA291619-rawdata.txt")

# remove RCM samples
df = df[,1:6]

# update colnames
colnames(df) = PRJNA291619_target$Sample

#create DGE class object
group <- as.factor(PRJNA291619_target$HeartFailure)
dge<- DGEList(counts=df, group=group)

#filter low expressed gene. Filtering keeps genes that have count-per-million (CPM) above k in n samples, where k is determined by min.count and by the sample library sizes and n is determined by the design matrix.
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization TMM. TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples.
dge <- calcNormFactors(dge)

# use limma voom
v <- voom(dge, plot=TRUE)

#save the R object
PRJNA291619_count = v$E
save(PRJNA291619_count, file = "HGEX_data/Study_Processing/PRJNA291619_count.ro")



