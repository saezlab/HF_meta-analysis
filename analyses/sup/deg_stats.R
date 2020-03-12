# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare the gene level statistics
#' of all different studies.

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

# 0. Processing object
experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

# p-values of each study
p_matrix = data.frame(get_all_limma(meta_list = METAheart,
                         limma_column = "P.Value"),
                      check.names = F) %>%
           mutate_all(function(x) -log10(x)) %>% 
           gather("experiment","pvalue") %>% na.omit()

# t-values of each study
t_matrix = data.frame(get_all_limma(meta_list = METAheart,
                         limma_column = "t"),
                      check.names = F) %>% 
  gather("experiment","tvalue") %>% na.omit()

# log fold change of each study
lfc_matrix = data.frame(get_all_limma(meta_list = METAheart,
                            limma_column = "logFC"),
                        check.names = F) %>% 
  gather("experiment","lfc") %>% na.omit()

deg_stats = list("pvalue" = p_matrix,
                 "tvalue" = t_matrix,
                 "lfc" = lfc_matrix)

saveRDS(deg_stats,
        "data/figure_objects/deg_stats.rds")





































































