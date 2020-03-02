# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Raw processing of external data set GSE84796

GSE84796_targets = t(read.table("data_processing/raw/GSE84796/GSE84796_rawtargets.txt", 
                                sep="\t",stringsAsFactors = F))

GSE84796_targets = GSE84796_targets[-1,c(1,9)] 
colnames(GSE84796_targets) = c("Sample","Disease")
GSE84796_targets = as_tibble(GSE84796_targets) %>% 
                   mutate("HeartFailure" = ifelse(grepl("end-stage",Disease),"yes","no"))


GSE84796_matrix = read.table("data_processing/raw/GSE84796/GSE84796_rawcounts.txt", 
                             sep = "\t", stringsAsFactors = F,
                             row.names = 1,header = T)

Ann_matrix = read.table("data_processing/raw/GSE84796/GPL15931-16833_red.txt", 
                        sep = "\t", stringsAsFactors = F,
                        header = T) %>% filter(GENE_SYMBOL!="")

Ann_matrix = filter(Ann_matrix, ID %in% rownames(GSE84796_matrix))
GSE84796_matrix = GSE84796_matrix[Ann_matrix$ID,GSE84796_targets$Sample]

#### Complete preproc ####
## mean by same probeset
GENENAMES = Ann_matrix$GENE_SYMBOL
HF_dataframe = GSE84796_matrix %>% mutate(ID = GENENAMES)
HF_dataframe = aggregate(x = HF_dataframe[, 1:(ncol(HF_dataframe))-1],
                         by = list(ID = HF_dataframe$ID), FUN = "mean", na.rm = T)
HF_expr_mat = as.matrix(HF_dataframe[,2:ncol(HF_dataframe)])
rownames(HF_expr_mat) = HF_dataframe$ID

HF_expr_mat = HF_expr_mat[,GSE84796_targets$Sample]

GSE84796_matrix = HF_expr_mat

save(GSE84796_matrix, file = "data_processing/processed/GSE84796_counts.ro")
save(GSE84796_targets, file = "data_processing/processed/GSE84796_targets.ro")



