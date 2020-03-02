# Author: Jan Lanzer, 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing Pipeline for Study Translatome (Spurrell et al 2019)
# License: MIT 

library(Biobase)
library(tidyverse)
library(limma)
library(edgeR)

# 1)  Read sample information, process, and save as target file. 
#read metainfo
translatometa= read.csv(file = "HGEX_data/Study_Processing/Metainformation/mmc1.csv", 
                        header = TRUE,
                        stringsAsFactors = F)

Translatome_targets = na.omit(translatometa[-c(81,82),c(2,4,6,7,8,9)])
colnames(Translatome_targets)[c(1,2,3,4,5,6)] = c("Sample","HeartFailure","Gender","SampleType",
                                              "Age","TechnicalTime")

Translatome_targets$HeartFailure = gsub("Control", "no", Translatome_targets$HeartFailure)
Translatome_targets$HeartFailure = gsub("DCM", "yes", Translatome_targets$HeartFailure)


Translatome_targets = as_tibble(apply(Translatome_targets,2,as.character))
Translatome_targets[78,"Gender"] = "f"

Translatome_targets = Translatome_targets %>% mutate(Gender = ifelse(Gender=="f",
                                                                     "female", "male"),
                                                     Age = as.numeric(Age))

colnames(Translatome_count) = gsub("_mR","",colnames(Translatome_count))
Translatome_count = Translatome_count[,Translatome_targets$Sample]

Translatome_targets$HTx = ifelse(grepl("LVAD",Translatome_targets$SampleType) & 
                                   !grepl("HTx",Translatome_targets$SampleType),
                                 "no","yes")

Translatome_targets$DCM = "yes"

technical_df = tibble(batch_reported = unique(Translatome_targets$TechnicalTime)) %>%
  mutate(batch_sign = LETTERS[1:length(batch_reported)])

Translatome_targets = left_join(Translatome_targets, technical_df,
                                by = c("TechnicalTime" = "batch_reported"))

## save meta
save(Translatome_targets , file = "HGEX_data/Study_Processing/Translatome_target.ro")



# 2) Read gene expression data, filter & normalize, and save as count file. 
# read GEX table
df = read.csv(file = "HGEX_data/Study_Processing/GEXraw/cardiacTranslatome_counts2.csv",
                        header = T,
                        stringsAsFactors = F)

#prefilter protein coding genes

df = df %>%
  filter(gene_biotype== "protein_coding") %>% 
  select(-gene_biotype) 

# Gene names should be in HGNC format, but contain some duplicates in this file (same Gene name but
# different ENS ID)
# Before selecting HGNC format, the mean of the expression value of those genes is calculated,
# and merged into one row.
  
dfmerge = df %>%
  group_by(gene_name) %>%
  select(-gene_id) %>%
  summarise_all(mean) %>% 
  ungroup() %>% 
  as.data.frame() %>%
  column_to_rownames("gene_name")

#create DGElistobject (no groups)
dge<- DGEList(counts=dfmerge)

#filter low expressed gene. Filtering keeps genes that have count-per-million (CPM) above k in n samples, where k is determined by min.count and by the sample library sizes and n is determined by the design matrix.
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization TMM. TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples.
dge <- calcNormFactors(dge)

# use limma voom
v <- voom(dge, plot=TRUE)

#save file 
Translatome_count = v$E

colnames(Translatome_count) = gsub("_mR","",colnames(Translatome_count))
Translatome_count = Translatome_count[,Translatome_targets$Sample]

save(Translatome_count, file = "HGEX_data/Study_Processing/Translatome_count.ro")

