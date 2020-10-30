# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script performs the functional analysis
# of the meta-ranked gene list, results are the ones
# provided in main figure

library(fgsea)
library(viper)
library(progeny)
library(WriteXLS)

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object
genesets = readRDS("data/Genesets_Dec19.rds") #MSigDB gene sets
# Load dororthea regulons
dorothea_regulon_human = 
  read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")
# We obtain the regulons based on interactions with confidence level A, B, C and D
regulons = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C","D")) %>%
  df2regulon()

# Iterators with experiments as names
experiments = names(METAheart)
names(experiments) = experiments

#a gene can be missing  in at most 4 experiments
fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10)
genes = names(fisher_rank)

lfc_matrix = (get_all_limma(meta_list = METAheart,
                            limma_column = "logFC"))[genes,]

# Generating ranks for GSEA
progeny_rank = dorothea_rank = gsea_rank = sort(rowMeans(sign(lfc_matrix),na.rm = T) * -log10(fisher_rank),
                 decreasing = T) #Ranking penalizing inconsistency in direction

saveRDS(progeny_rank, 
        file = "data/shiny/directed_ranking.rds")

#1. GSEA
red_gensets = genesets[c("MSIGDB_CANONICAL","MSIGDB_HMARKS",
                               "MSIGDB_GO_CELLCOMP","MSIGDB_GO_BIOLPROC",
                               "MSIGDB_GO_MOLFUNC")]

gsea_meta = lapply(red_gensets, function(x){
  set.seed(1234) # fgsea unstable results
  GSEA_results = fgsea(pathways = x, stats = gsea_rank,
                       minSize = 15, maxSize = 300, 
                       nperm = 1000) %>% as_tibble() %>% 
    arrange(desc(abs(NES)))
})

#For plotting
selected_ps = enframe(gsea_meta) %>% unnest() %>% 
  arrange(desc(abs(NES))) %>% dplyr::slice(1:50) %>% filter(!duplicated(pathway)) %>%
  arrange(NES) %>% mutate(pathway = strtrim(pathway,40))

selected_ps = selected_ps[!duplicated(selected_ps$pathway),]

saveRDS(selected_ps, 
        file = "data/figure_objects/GSEA_results_sel.rds")

#For supplementary
all_GSEA = enframe(gsea_meta,name = "database","value") %>% unnest() %>% 
  arrange(desc(abs(NES))) %>% dplyr::mutate(padj = p.adjust(pval,method = "BH"))

print("N gene sets")
print(dim(all_GSEA))

print("N p_val<0.05")
print(dim(all_GSEA %>% dplyr::filter(pval<=0.05)))

saveRDS(all_GSEA, 
        file = "data/shiny/GSEA_results.rds")

# 1.1 miRNAs
gset = genesets$MSIGDB_MIRNA
set.seed(1234) # fgsea unstable results

#miRNA_results = fgsea(pathways = gset, stats = gsea_rank,
#                     minSize = 15, maxSize = 300, nperm = 1000) %>% as_tibble() %>% 
#  arrange(desc(abs(NES))) %>% dplyr::select(pathway, pval, padj, ES, NES)

# What if we try viper

miRNAs_regulons = enframe(gset,value = "target",name = "tf") %>% 
  unnest() %>% mutate(mor = 1, likelihood = 1) %>%
  df2regulon()

miRNA_results = msviper_summary(msviper(gsea_rank,
                                              miRNAs_regulons,
                                              minsize = 10,
                                              verbose = FALSE))

print("N miRNA")
print(dim(miRNA_results))

print("N p_val<0.05")
print(dim(miRNA_results %>% dplyr::filter(pvalue<=0.05)))

saveRDS(miRNA_results, 
        file = "data/figure_objects/GSEA_mir_results.rds")

saveRDS(miRNA_results, 
        file = "data/shiny/GSEA_mir_results.rds")

#2. Dorothea
dorothea_results = msviper_summary(msviper(dorothea_rank,
                                           regulons,
                                           minsize = 10,
                                           verbose = FALSE))

print("N TFs")
print(dim(dorothea_results))

print("N p_val<0.05")
print(dim(dorothea_results %>% dplyr::filter(pvalue<=0.05)))

saveRDS(dorothea_results, 
        file = "data/figure_objects/dorothea_results.rds")

saveRDS(dorothea_results, 
        file = "data/shiny/dorothea_results.rds")

#3. PROGENy
progeny_rank_mat = as.matrix(progeny_rank)
colnames(progeny_rank_mat) = "meta"

meta_progeny = t(progeny(progeny_rank_mat,scale = FALSE, organism = "Human",
                         top = 200, perm =1))

#Calculation of pvalues permuting rank
set.seed(1234)

permutation_prog = sapply(1:1000, function(x){
  
  rnd_order = progeny_rank_mat
  rownames(rnd_order) = sample(rownames(rnd_order))
  progeny(rnd_order,scale = FALSE, organism = "Human",
          top = 200, perm =1)
  
})

rownames(permutation_prog) = rownames(meta_progeny)

progeny_pvals = meta_progeny[,1]

for(p in names(meta_progeny[,1])){
  
  metascore = meta_progeny[p,1]
  
  if(sign(metascore)>0){
    
    pval = sum(permutation_prog >= metascore)/length(permutation_prog)
    
  }else{
    
    pval = sum(permutation_prog <= metascore)/length(permutation_prog)
    
  }
  
  progeny_pvals[p] = pval
  
}

prog_res = tibble(progeny_scores = meta_progeny[,1], 
                  progeny_pvals, pathway = names(progeny_pvals)) %>%
           arrange(progeny_pvals)

print("N paths")
print(dim(prog_res))

print("N p_val<0.05")
print(dim(prog_res %>% dplyr::filter(progeny_pvals<=0.05)))


saveRDS(prog_res, 
        file = "data/figure_objects/PROGENy_results.rds")

saveRDS(prog_res, 
        file = "data/shiny/PROGENy_results.rds")

#Generating supplementary table


WriteXLS(x = c("all_GSEA",
               "dorothea_results",
               "prog_res",
               "miRNA_results"), 
         ExcelFileName = "data/paper_sup/SupplementalTable3.xlsx",
         SheetNames = c("GSEA",
                        "TF_activities",
                        "Pathway_activities",
                        "miRNAs")
)


























