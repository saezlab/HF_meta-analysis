# Author: Jan Lanzer 2019
# Contact: jan.lanzer@bioquant.uni-heidelberg.de
# Description: Processing pipeline for Study GSE123976
# License: MIT 

library(tidyverse)  
library(limma)
library(edgeR)
library(stringr)


# 1 ) Read sample information, process, and save as target file. 
# load dictionary between SRR and GSM sample IDs
dictionary =read.table("HGEX_data/Study_Processing/Metainformation/PRJNA510407.txt", 
                       header = T, 
                       stringsAsFactors = F) %>%
  rename(GSM = submitted_ftp)

# load meta information from GEO
meta = read.csv("HGEX_data/Study_Processing/Metainformation/GSE123976_information.txt_edit.txt",
                sep = "\t",
                skip = 3,
                header = TRUE,
                stringsAsFactors = F) %>%
  as.tibble %>%
  select(Title, Characteristics, Library.source, GSM) %>% 
  inner_join(.,dictionary) %>% 
  filter(Library.source == "transcriptomic") # here we filter out bisulfite seq
  
# extract age and gender information
feat = gsub(".* ","", str_split_fixed(meta$Characteristics, ",", 6)) 

GSE123976_target = meta %>% full_join(data.frame(feat, "GSM"= meta$GSM))

colnames(GSE123976_target) = c(colnames(GSE123976_target[1:5]), "Biopsy", "HeartFailure", "Age", "Gender", "Ethnicity", "RNA")

GSE123976_target = GSE123976_target%>%
  as.tibble() %>%
  select(run_accession, HeartFailure, Age, Gender)  %>%
  mutate(HeartFailure = as.character(HeartFailure),
         Age = as.numeric(as.character(Age)), 
         Gender = as.character(Gender), 
         DCM = c(rep("no",3),"no", "yes", "yes", "yes", "no", "no"), # this information is taken from the publication https://journals.physiology.org/doi/full/10.1152/ajpheart.00016.2019
         HTx = "no")  %>%  # all samples were acquired during LVAD  
  rename(Sample = run_accession)

GSE123976_target$HeartFailure = gsub("CON", "no", GSE123976_target$HeartFailure)
GSE123976_target$HeartFailure = gsub("HF", "yes", GSE123976_target$HeartFailure)

GSE123976_target$Gender = gsub("F", "female", GSE123976_target$Gender)
GSE123976_target$Gender = gsub("M", "male", GSE123976_target$Gender)

# save r-object
save(GSE123976_target, file = "HGEX_data/Study_Processing/GSE123976_target.ro")


#2 ) Read gene expression data, filter & normalize, and save as count file. 
# read GEX data table (as downloaded from biojupies)
df= read.table("HGEX_data/Study_Processing/GEXraw/GSE123976-rawdata.txt")

#check if sample names are in the same order as in the target file 
colnames(df) == GSE123976_target$Sample

#create DGE class object
group <- GSE123976_target$HeartFailure
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization method TMM (trimmed mean of M)
dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

GSE123976_count= v$E
save(GSE123976_count, file = "HGEX_data/Study_Processing/GSE123976_count.ro")

