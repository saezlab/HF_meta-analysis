# Author: Jan Lanzer, 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing Pipeline for Study PRJNA477855
# License: MIT 

library(Biobase)
library(tidyverse)
library(limma)
library(edgeR)

# 1)  Read sample information, process, and save as target file. 
# read meta table
meta = read.table("HGEX_data/Study_Processing/Metainformation/PRJNA477855_information_edit.txt",sep = "\t", skip = 3, header = TRUE)

# select relevant columns in meta table
meta = meta %>% as_tibble() %>% select(GSM, Title, Source.name, Characteristics, Platform.ID)

# split the "Characteristics column into single pieces of information
info2= t(matrix(unlist(strsplit(gsub(" ","", meta$Characteristics),",")),nrow = 4))

# remove unecessary infos
info2[,4]=gsub("Sex:","", info2[,4])
info2[,3]=gsub("age:","", info2[,3])
info2[,2]=gsub("disease:","", info2[,2])

info2 = info2 %>%as.tibble() %>% mutate(Sample = meta$Title) %>% dplyr::rename(Disease = V2)

# read meta info from paper supplement
cmeta = as.tibble(read.csv(file = "HGEX_data/Study_Processing/Metainformation/PRJNA477855_clinicalmeta.csv", skip = 1))

# merge all three metainformation tables, select relevant and non redundant columns and format values
target = info2 %>% 
  inner_join(cmeta, by = "Sample") %>% 
  select(-V1, -V3, -V4) %>% 
  inner_join(meta %>% 
  rename(Sample = Title) %>% 
  select(GSM,Sample)) %>%
  rename(HeartFailure = Sample) %>% 
  rename(Sample = GSM) %>% 
  rename(Age = Age.at.Transplant) %>% 
  rename(Gender = Sex) %>%
  mutate(HeartFailure = "yes") %>% 
  mutate(Gender = tolower(Gender)) %>%
  mutate(Sample = as.character(Sample))


target$HeartFailure[str_detect(target$Disease, "non-failing")] = "no"
target = target %>% select(-Race, -Ethnicity,-Cause.of.Death..NF.only.,-RNA.integrity..RIN.) %>%
  mutate(DCM = "no") %>% 
  mutate(HTx = "yes")

target$DCM[str_detect(target$Disease, "dilated")] = "yes"

PRJNA477855_target = target  %>% select(-Disease)


# save Targetfile and GEXtable (normalized, filtered, in logCPMM)
save(PRJNA477855_target, file = "HGEX_data/Study_Processing/PRJNA477855_target.ro")



# 2) Read gene expression data, filter & normalize, and save as count file. 
# read GEX table
df = read.table("HGEX_data/Study_Processing/GEXraw/PRJNA477855-rawdata.txt")
colnames(df) = PRJNA477855_target$Sample

#create DGE class object
group <- as.factor(PRJNA477855_target$HeartFailure)
dge<- DGEList(counts=df, group=group)

#filter low expressed gene. Filtering keeps genes that have count-per-million (CPM) above k in n samples, where k is determined by min.count and by the sample library sizes and n is determined by the design matrix.
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization TMM. TMM is the recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples.
dge <- calcNormFactors(dge)

# use limma voom
v <- voom(dge, plot=TRUE)

#save the R object
PRJNA47855_count = v$E
save(PRJNA47855_count, file = "HGEX_data/Study_Processing/PRJNA477855_count.ro")
