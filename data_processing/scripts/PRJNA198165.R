# Author: Jan Lanzer 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing pipeline for Study PRJNA198165
# License: MIT 

library(tidyverse)  
library(limma)
library(edgeR)
library(Biobase)
library(plyr)


# 1 ) Read sample information, process, and save as target file. 
#### read target information
metaNF = read.csv("HGEX_data/Study_Processing/Metainformation/PRJNA198165_clinicalmeta.csv", 
                  header = TRUE, 
                  skip = 1, 
                  stringsAsFactors = F)
metaCM = read.csv("HGEX_data/Study_Processing/Metainformation/NIHMS559217-supplement-Supplemental_Material-24-25-1_mod.csv", 
                  header = TRUE, 
                  skip = 0,
                  stringsAsFactors = F)

# rename column titles and merge both dataframes
metaNF = metaNF[,1:7] %>%
  dplyr::rename(Title = ID) %>% 
  dplyr::rename(History.of.Smoking = Smoker..yr.)

metaCM= metaCM %>% 
  as_tibble() %>% 
  dplyr::rename(Age= Age..yr.)

PRJNA198165_target = rbind.fill(metaNF, metaCM)

# delete postLVAD_Samples, they are replicates from the same patients as the preLVAD_Samples and 
# shall not be used in this meta analysis
PRJNA198165_target= PRJNA198165_target[1:24,]

##### create target file
group <- c(rep(1,8),rep(2,16))

# Structure and rename information according to our convention 
PRJNA198165_target = PRJNA198165_target %>% 
  dplyr::rename(Sample = Title) %>%
  dplyr::rename(Gender = Sex) %>% 
  mutate(HeartFailure = group) %>%
  as.tibble()

PRJNA198165_target$Gender= gsub("M", "male", PRJNA198165_target$Gender)
PRJNA198165_target$Gender= gsub("F", "female", PRJNA198165_target$Gender)
PRJNA198165_target$HeartFailure= gsub(1, "no", PRJNA198165_target$HeartFailure)
PRJNA198165_target$HeartFailure= gsub(2, "yes", PRJNA198165_target$HeartFailure)

PRJNA198165_target = PRJNA198165_target %>% 
  select(Sample, HeartFailure, Age, Gender, Diagnosis) %>% 
  mutate(HTx = "no") %>%
  mutate(DCM = "no")

PRJNA198165_target[PRJNA198165_target$HeartFailure == "no", "Diagnosis"] = "control"
PRJNA198165_target[PRJNA198165_target$Diagnosis == "NICM", "DCM"] = "yes"
PRJNA198165_target[PRJNA198165_target$HeartFailure =="no", "HTx"] = "yes"

PRJNA198165_target = PRJNA198165_target %>% select(-Diagnosis)

#save targetfile
  save(PRJNA198165_target, file = "HGEX_data/Study_Processing/PRJNA198165_target.ro")



#### 2) Read gene expression data, filter & normalize, and save as count file. 

#### read GEX data table (as downloaded from biojupies)
df = read.table("HGEX_data/Study_Processing/GEXraw/PRJNA198165-rawdata.txt")

# reordering samples, sorting samples by their #ID
colnames(df)= gsub("sample", "", colnames(df))
df= df[,order(strtoi(colnames(df)))]
colnames(df)= PRJNA198165_target$Sample

# delete postLVAD_Samples, they are replicates from the same patients as the preLVAD_Samples and shall not be used in this meta analysis
df=df[,1:24]

#### Filter & Normalize
#create DGE class object
dge<- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

PRJNA198165_count= v$E
save(PRJNA198165_count, file = "HGEX_data/Study_Processing/PRJNA198165_count.ro")
