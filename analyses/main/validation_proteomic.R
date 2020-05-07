# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jan D. Lanzer]
# roramirezf@uni-heidelberg.de

#' Description: Here we trace the cardiac origin
#' of HF plasma proteomic biomarkers from Egerstedt (PMID: 31862877)
#' 
#' 

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
library(biomaRt)
library(tidyr)
library(tibble)
library(dplyr)
library(fgsea)

## prepare conensus signature 
METAheart = readRDS(file = "data/METAheart.rds") #main object
meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample,HeartFailure,ExpID) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_")) 

fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10)

genes = names(fisher_rank)
  
# p-values of each study
p_matrix =(get_all_limma(meta_list = METAheart,
                         limma_column = "P.Value"))[genes,]

# t-values of each study
t_matrix =(get_all_limma(meta_list = METAheart,
                         limma_column = "t"))[genes,]
  
# log fold change of each study
lfc_matrix = (get_all_limma(meta_list = METAheart,
                            limma_column = "logFC"))[genes,]

#### Analyze Proteome 
processPROTEOMICS = function(filepath = "data/proteomics/feature_dict.tsv") {
  con = file(filepath, "r")
  line_df = c()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    split_line = unlist(strsplit(line,split="\t"))
    split_line = sapply(split_line,function(x) gsub("\"","",x))
    line_df = rbind(line_df,split_line)
  }
  close(con)
  return(line_df)
}

proteins_dictionary = data.frame(processPROTEOMICS()[-1,c(2,3)],
                                 stringsAsFactors = F) %>% as_tibble() 

colnames(proteins_dictionary) = c("protein","uniprot")

proteins_dictionary$uniprot = sapply(proteins_dictionary$uniprot, 
                                     function(x) unlist(strsplit(x,", ")))

proteins_dictionary = proteins_dictionary %>% unnest()

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

g_annotation = getBM(filters= "uniprotswissprot", 
                     attributes= c("description","uniprotswissprot","hgnc_symbol"),
                     values=proteins_dictionary$uniprot,mart= mart,uniqueRows = F)

proteins_dictionary = left_join(proteins_dictionary, g_annotation, 
                                by = c("uniprot" = "uniprotswissprot"))

proteins_dictionary = na.omit(unique(proteins_dictionary))



#### Reading result tables

# Manifest HF 
manifest_res = as_tibble(data.frame(processPROTEOMICS("data/proteomics/manifestHF.tsv")[-1,1:5],
                                    stringsAsFactors = F))

colnames(manifest_res) = c("protein","or","beta","se","pval")

manifest_res = manifest_res %>%
  mutate_at(c("or","beta","se","pval"), as.numeric) %>%
  left_join(proteins_dictionary) %>% 
  dplyr::mutate(adj_pval = p.adjust(p = pval,method = "BH")) %>%
  na.omit()

plasma_BM_manifest = manifest_res %>% dplyr::filter(adj_pval<0.01,
                                                    hgnc_symbol %in% genes)

plasma_BM_manifest$meta_rank = unlist(sapply(plasma_BM_manifest$hgnc_symbol, 
                                             function(x) which(genes %in% x)))

# Development HF 
dev_res = as_tibble(data.frame(processPROTEOMICS("data/proteomics/developmentHF.tsv")[-1,1:5],
                                    stringsAsFactors = F))

colnames(dev_res) = c("protein","hr","beta","se","pval")

dev_res = dev_res %>%
  mutate_at(c("hr","beta","se","pval"), as.numeric) %>%
  left_join(proteins_dictionary) %>% 
  dplyr::mutate(adj_pval = p.adjust(p = pval,method = "BH")) %>%
  na.omit()

plasma_BM_dev = dev_res %>% dplyr::filter(adj_pval<0.01, 
                                          hgnc_symbol %in% genes)

#save a copy without pval
plasma_BM_dev_full = dev_res %>% dplyr::filter(hgnc_symbol %in% genes)

plasma_BM_dev$meta_rank = unlist(sapply(plasma_BM_dev$hgnc_symbol, 
                                        function(x) which(genes %in% x)))

# save Biomarker list for furhter analysis in validation_plotting.R
saveRDS(list("manifest" = unique(plasma_BM_manifest$hgnc_symbol),
     "development" = unique(plasma_BM_dev$hgnc_symbol)), "data/figure_objects/protein_genesets.rds")
 


## Inference of Biomarkers of myocardial origin - MANIFEST HF
meant = rowMeans(t_matrix,
                 na.rm = T)

manifest_agrmnt = plasma_BM_manifest %>% 
  dplyr::filter(hgnc_symbol %in% names(meant),
                meta_rank <= 3000) %>%
  dplyr::select(hgnc_symbol,or) %>% 
  dplyr::mutate(logor = log2(or)) %>% 
  unique()

manifest_agrmnt$trnscrpt_evid_t = meant[manifest_agrmnt$hgnc_symbol]

# save Biomarker list for plotting in validation_plotting.R
saveRDS(manifest_agrmnt, "data/figure_objects/proteomics_manifest.rds")



## Inference of Biomarkers of myocardial origin - DEVELOPMENT HF
meant = rowMeans(t_matrix,
                 na.rm = T)

dev_agrmnt = plasma_BM_dev_full %>% 
  dplyr::filter(hgnc_symbol %in% names(meant))  %>%
  dplyr::select(hgnc_symbol,hr) %>% 
  dplyr::mutate(loghr = log2(hr)) %>% 
  unique()

dev_agrmnt$trnscrpt_evid_t = meant[dev_agrmnt$hgnc_symbol]

# save Biomarker list for plotting in validation_plotting.R
saveRDS(dev_agrmnt, "data/figure_objects/proteomics_dev.rds")

