# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script makes the meta-analysis of all
# studies and measures the contribution of each one in the
# final ranking. Correlates contribution with sample size
# and 

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library('biomaRt')
library(WriteXLS)
library(fgsea)

# 0. Processing object

experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

# Generating object for barplots

meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample,HeartFailure,ExpID) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_")) 

saveRDS(meta_targets, 
        file = "data/figure_objects/meta_targets.rds")


# 1. Meta-analysis

fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10)

print("Total number o")
print(length(fisher_rank))
sum(fisher_rank < .00005)

genes = names(fisher_rank)

shiny_allgenes = get_tibble_union(meta_list = METAheart,
                                  index_name = "HF_limma") %>% 
                 filter(ID %in% genes)

write.table(shiny_allgenes,col.names = T,
            row.names = F,sep = "\t",quote = F,
            file = "data/shiny/all_summary_stats.txt")

saveRDS(fisher_rank, 
        file = "data/shiny/fisher_rank.rds")

# p-values of each study
p_matrix =(get_all_limma(meta_list = METAheart,
                         limma_column = "P.Value"))[genes,]

# t-values of each study
t_matrix =(get_all_limma(meta_list = METAheart,
                         limma_column = "t"))[genes,]

# log fold change of each study
lfc_matrix = (get_all_limma(meta_list = METAheart,
                            limma_column = "logFC"))[genes,]

# 1.2 print, where are the markers

marker_genes = c("MYH6","MME","CNN1","NPPA","KCNH2","SLC2A1",
                 "ATP2A2","COL21A1","COL15A1","ECM2","MXRA5",
                 "KIT","FNDC1","LAMA4","SSPN","KCNN3","FGF14")

print(data.frame(Rank = which(genes %in% marker_genes),
           mgens = genes[which(genes %in% marker_genes)]))

# 2. Calculate the contribution of each study

meta_ranking = -log10(fisher_rank)

study_deg_list = lapply(METAheart, function(x){
  deg = dplyr::slice(x$HF_limma, #top 500 genes of each study
                     1:500)
  return(deg[[1]])
})

set.seed(1234)
contribution = fgsea(pathways = study_deg_list,
                     stats = meta_ranking,nperm = 1000) %>% 
  as_tibble()

print(contribution,n=50)

saveRDS(contribution, file = "data/figure_objects/contribution.rds")

# 3. No correlation between contribution and sample size 

cor_contribt_size = cor.test(experiment_size[contribution$pathway],
                             contribution$ES, method = "spearman")

print("correlation: contribution-size")

print(cor_contribt_size)

# 2. Annotate genes for summary
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

g_annotation = getBM(filters= "hgnc_symbol", 
                     attributes= c("description","hgnc_symbol"),
                     values=genes,mart= mart)

# 3. Generate summary object 
fisher_summary = tibble("gene" = genes,
                        "fisher_pvalue" = fisher_rank,
                        "mean_t" = rowMeans(t_matrix,na.rm = T),
                        "mean_lfc" = rowMeans(lfc_matrix, na.rm = T)) %>% 
  left_join(g_annotation, by = c("gene"="hgnc_symbol")) %>%
  dplyr::select(gene, description, fisher_pvalue,mean_t,
                mean_lfc)

write.table(fisher_summary,col.names = T,
            row.names = F,quote = F,
            sep = "\t",
            file = "data/meta_analysis_summary.txt")

p_matrix = as.data.frame(p_matrix) %>% rownames_to_column("gene")
t_matrix = as.data.frame(t_matrix) %>% rownames_to_column("gene")
lfc_matrix = as.data.frame(lfc_matrix) %>% rownames_to_column("gene")

WriteXLS(x = c("fisher_summary",
               "p_matrix",
               "t_matrix",
               "lfc_matrix"), 
         ExcelFileName = "data/paper_sup/meta_analysis.xlsx",
         SheetNames = c("meta_analysis_summary",
                        "individual_pvalues",
                        "individual_tvalues",
                        "individual_lfc")
)






