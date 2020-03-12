# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the robustness of
#' our replicability measurements

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(tidyverse)
library(cowplot)

# For labeling
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)


# 1. Disease score

pairwise_auc = readRDS(file = "data/figure_objects/robust_ds.rds")

pairwise_auc = filter(pairwise_auc, PredictorExperiment!=PredictedExperiment)

pairwise_plot = ggplot(pairwise_auc, aes(x = factor(n_genes,
                                                    levels = c("50","100","200",
                                                               "500","1000")),
                                         y = single,
                                         color = PredictorExperiment)) +
  stat_boxplot() + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=1)) +
  xlab("Number of genes in classifier") +
  ylab("AUC of Disease Score") + labs(color = "Predictor")

# 2. Enrichment score

pairwise_es_res= readRDS(file =  "data/figure_objects/robust_es.rds")

pairwise_es_res = filter(pairwise_es_res, 
                         Reference!=DEG) %>%
                  mutate(direction = factor(direction,
                                            levels = c("up","down")))

pairwise_plot_es = ggplot(pairwise_es_res, aes(factor(n_genes,
                                                      levels = c("50","100","200",
                                                                 "500","1000")),
                                         y = ES,
                                         color = factor(DEG,
                                                        levels = names(experiment_size)))) +
  stat_boxplot() + 
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=1)) +
  xlab("Number of DEG chosen") +
  ylab("Enrichment Score in background experiments") +
  facet_grid(direction ~.) + labs(color = "DEG")


# Joint figure
pdf("./analyses/figures/sup/robustness_es_ds.pdf",
    width = 16,
    height = 8)

plot(plot_grid(pairwise_plot,pairwise_plot_es, 
          nrow = 1,rel_widths = c(1,1),labels = c("A","B")))

dev.off()













