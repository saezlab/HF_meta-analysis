# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script test the best number of genes from the meta-analysis
# using the classification performance of disease score using different gene lists

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

# 0. Processing object

experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

# 1. Meta-analysis

fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10)

t_matrix = get_all_limma(meta_list = METAheart,
                         limma_column = "t")
lfc_matrix = get_all_limma(meta_list = METAheart,
                           limma_column = "logFC")
genes = names(fisher_rank)

# 2. How the performance improves as we include more and more top genes

ngenes = seq(50,14000,100)
names(ngenes) = ngenes

InPerformance = map(ngenes, function(ingenes){ # We will include 100 genes at a time til 5k
  testgenes = names(fisher_rank[1:ingenes]) # Defining the genes to be used to calculate DS
  ingenes_stats = getRisk_Stats_v2(Experiment_List = METAheart,
                                   limma_t_mat = t_matrix, 
                                   genes = testgenes)
  
  ingenes_stats = map(ingenes_stats, function(x){ #Calculating mean AUC of each predicted
    enframe(x[["SingleAUC"]],       #experiment
            "PredictorExperiment",
            "AUC")
  }) %>% enframe("PredictedExperiment") %>%
    unnest() 
  #%>% group_by(PredictedExperiment) %>%
  #summarise(meanAUC = mean(AUC))
  
  return(ingenes_stats)
}) %>% enframe("nin") %>% unnest() %>% mutate(nin = as.numeric(nin),
                                              PredictedExperiment = factor(PredictedExperiment,
                                                                           levels = names(experiment_size)
                                              ))

InPerformance_df = InPerformance %>% group_by(PredictedExperiment,nin) %>%
  summarise(meanAUC = mean(AUC))

# Best performance: Saturation

saveRDS(InPerformance_df, 
        file = "data/figure_objects/InPerformance_df.rds")


InPerformance_zoom_df = InPerformance_df

max_class = InPerformance_zoom_df %>% 
  arrange(desc(meanAUC)) %>% 
  slice(1)

sum(max_class$nin <= 500)

saveRDS(InPerformance_zoom_df, 
        file = "data/figure_objects/InPerformance_df_zoom.rds")

# 3. How the performance worsens as we exclude top genes

ngenes = seq(50,14000,100)
names(ngenes) = ngenes

OutPerformance = map(ngenes, function(outgenes){ # We will include 100 genes at a time til 5k
  testgenes = names(fisher_rank[outgenes:length(fisher_rank)]) # Defining the genes to be used to calculate DS
  ingenes_stats = getRisk_Stats_v2(Experiment_List = METAheart,
                                   limma_t_mat = t_matrix, 
                                   genes = testgenes)
  
  ingenes_stats = map(ingenes_stats, function(x){ #Calculating mean AUC of each predicted
    enframe(x[["SingleAUC"]],       #experiment
            "PredictorExperiment",
            "AUC")
  }) %>% enframe("PredictedExperiment") %>%
    unnest() %>% group_by(PredictedExperiment) %>%
    summarise(meanAUC = mean(AUC))
  
  return(ingenes_stats)
}) %>% enframe("nout") %>% unnest() %>% 
  mutate(nout = as.numeric(nout),
         PredictedExperiment = factor(PredictedExperiment,
                                      levels = names(experiment_size)))

saveRDS(OutPerformance, 
        file = "data/figure_objects/OutPerformance_df.rds")

