# Author: Jan Lanzer 2020
# Description: This script performs an analysis for fetal samples from the meta heart and compare
#               Results with those obtained from the meta heart analysis
  # fetal study: GSE52061, PRJNA522417

library(tidyverse)
library(limma)
library(fgsea)
library(viper)
library(progeny)
library(WriteXLS)
library("Hmisc")
library(ggrepel)

#files
fetal_experiments = readRDS("HGEX_data/fetal_METAheart.rds")
METAheart= readRDS("HGEX_data/METAheart.rds")
meta_rank = as_tibble(read.csv("HGEX_data/METArank_March2020.csv",
                               header = T,
                               sep= ",",
                               stringsAsFactors = F) )
#source code
source(file = "HGEX_src/data_utils.R")
source("HGEX_src//misc_utils.R")


### 2 fetal studies can be analyzed
#select which study should be analyzed

#study= "GSE52601"
study= "PRJNA522417"

targets = fetal_experiments[[study]]$TARGETS
count = fetal_experiments[[study]]$GEX

###

fetal = as.factor(targets$fetal)
design = model.matrix(~fetal)
rownames(design) <- targets$Sample
design

fit = lmFit(count, design)
fit2 = eBayes(fit)

fetalDEA = as.data.frame(topTable(fit2,adjust.method = "BH",coef = "fetalyes",
                                  number = Inf)) %>%
  rownames_to_column("gene") %>%
  arrange(desc(abs(t))) %>% as_tibble()

# volcanoplot of DEA
ggplot(data= fetalDEA, mapping = aes(x=logFC, y= -log10(adj.P.Val)))+
  geom_point()+
  ggtitle("fetal vs adult non failing")

#save for further plotting in Script validation_plotting.R
saveRDS(fetalDEA, file =paste0("HGEX_data/figure_objects/fetalDEgenes_",study,".rds"))

### Compare DE genes from fetal study with top500 genes in meta ranking:

# volcanoplot of DEA and colored genes from the top500 in HF consensus signature
fetalDEA_mod = fetalDEA %>% 
  mutate(metaTop500= gene %in% meta_rank$gene[1:500]) %>% 
  arrange(metaTop500)
  


ggplot(data= fetalDEA_mod, mapping = aes(x=logFC, y= -log10(adj.P.Val) ))+
  geom_point(aes(color= metaTop500))+
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = 2)+
  ggtitle("fetal vs adult non failing")+
    theme_minimal()+
  geom_label_repel(aes(label=ifelse( (metaTop500 == T) & (-log10(adj.P.Val)> 5), gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') 




### Enrichment Analysis with fgsea
## geneset and ranking preparation
# calculate -log10 pvalue
meta_rank = meta_rank %>%
  mutate(log10pval = -log10(fisher_pvalue)) %>%
  arrange(desc(log10pval))

# filter DE genes from fetal samples for significance
fetalDEA_sig= fetalDEA %>%
  filter(adj.P.Val< 0.05) %>% arrange(adj.P.Val)

# create directional gene lists
fetalDEA_up = fetalDEA_sig %>% filter(logFC > 0)
fetalDEA_dn = fetalDEA_sig %>% filter(logFC < 0)


## 1. Enrichment
## Enrich all DE fetal genes in the -log10(fisher_p) sorted meta rank (undirected)
geneset = list(fetalDEA_sig$gene)
gsea_rank_undir= meta_rank$log10pval
names(gsea_rank_undir) = meta_rank$gene
fgseaRes_undir = fgsea(pathways=  geneset,
                       stats = gsea_rank_undir,
                       nperm= 1000,
                       minSize = 2,
                       gseaParam = 2) %>%
  as_tibble()

#save for further plotting in Script validation_plotting.R
saveRDS(geneset,file =  paste0("HGEX_data/figure_objects/fealgeneset",
                        study,
                        ".rds"))

## 2. Enrichment
# Enrich up and down regulated DE genes from fetal separatley in a directed meta ranking
# create directed metaranking
# lfc_matrix of top differentially expressed genes from the meta ranking in all studies
lfc_matrix = (get_all_limma(meta_list = METAheart,
                            limma_column = "logFC"))[meta_rank$gene,]

# create a gene level statistic that penalizes for inconsistency in direction
gsea_rank = sort(rowMeans(sign(lfc_matrix),na.rm = T) * -log2(meta_rank$fisher_pvalue),
                 decreasing = T) #Ranking penalizing inconsistency in direction

# up regulated and down regulated fetal genes are two gene sets
geneset= list("up" =fetalDEA_up$gene,"dn"= fetalDEA_dn$gene)

# enrichment performed
fgseaRes_dir = fgsea(pathways = geneset,
                     stats = gsea_rank,
                     nperm = 1000,
                     minSize = 2,
                     gseaParam = 2) %>%
  as_tibble()%>%
  arrange(desc(NES))


### Transcription factor analysis

##Function to group Dorothea regulons. (Can be deleted when updated source script is used)
## Input: A data frame containing Dorothea regulons, as stored in
## https://github.com/saezlab/ConservedFootprints/tree/master/data
## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}

## Load dororthea regulons
dorothea_regulon_human =
  read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")

# We obtain the regulons based on interactions with confidence level A, B and C
regulons = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C","D")) %>%
  df2regulon()

# TF activities from fetal
dorothea_rank = fetalDEA$t
names(dorothea_rank)= fetalDEA$gene

dorothea_results_fet = msviper_summary(msviper(dorothea_rank,
                                           regulons,
                                           minsize = 10,
                                           verbose = FALSE))

#TF activities from metaheart
dorothea_results_meta = msviper_summary(msviper(gsea_rank,
                                           regulons,
                                           minsize = 10,
                                           verbose = FALSE))


#compare results
setfetal = dorothea_results_fet %>% filter(pvalue <0.05)
setmeta = dorothea_results_meta %>% filter(pvalue<0.05)

#merged list of significant TFs
sig_TF = setfetal %>% 
  inner_join(setmeta, by = "RegulonName")
dim(sig_TF) #13 tf with p val < 0.05  ; 0 TFs with adj pvalue <0.05

## check the correlation 
#corrleation of full merged TF list
full_TF = dorothea_results_fet %>% 
  inner_join(dorothea_results_meta, by = "RegulonName")

fullTF_red = full_TF %>%
  dplyr::select(RegulonName,NES.x, NES.y) %>% column_to_rownames("RegulonName")

corFullTF = rcorr(as.matrix(fullTF_red))

#correlation of only significant TFs
TF_red =sig_TF %>% dplyr::select(RegulonName,NES.x, NES.y) %>% column_to_rownames("RegulonName")
corrSigTF = rcorr(as.matrix(TF_red))

#List of merged significant TFs with SAME direction:
targets = TF_red %>% 
  rownames_to_column("RegulonName") %>% 
  filter(sign(NES.x) == sign(NES.y))%>%
  mutate(corr = cor(NES.x,NES.y))

#plotting table
plot.TFs.data = dorothea_results_fet %>% 
  left_join(dorothea_results_meta, by = "RegulonName") %>% 
  filter(pvalue.x <0.05) %>% 
  mutate(quadrant = (sign(NES.x) == sign(NES.y)),
                         metasig = (pvalue.y <0.05),
                         colored = (quadrant == T) & (metasig == T)) %>%
  rename(NES.fetal = NES.x, 
         NES.meta = NES.y)

#save plotting table 
saveRDS(plot.TFs.data, file = paste0("HGEX_data/figure_objects/fetalTFs_",study, ".rds"))
  
