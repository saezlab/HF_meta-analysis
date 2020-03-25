# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the statistics
#' of limma from different studies

library(tidyverse)
library(cowplot)

METAheart = readRDS(file = "data/METAheart.rds") #main object
experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

deg_stats = readRDS("data/figure_objects/deg_stats.rds")

pval = ggplot(deg_stats$pvalue,
       aes(x = factor(experiment,levels = names(experiment_size)), 
           y = pvalue)) + geom_violin() +
  theme_minimal() +
  theme(axis.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=10),
        panel.background = element_rect(fill=NULL,
                                        color = "black"),
        panel.grid = element_blank()) +
        xlab("dataset") + ylab("-log10(p-value)")

tval = ggplot(deg_stats$tvalue,
       aes(x = factor(experiment,levels = names(experiment_size)), 
           y = tvalue)) + geom_violin() +
  theme_minimal() +
  theme(axis.title = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=10),
        panel.background = element_rect(fill=NULL,
                                        color = "black"),
        panel.grid = element_blank()) +
  xlab("dataset")

lfc = ggplot(deg_stats$lfc,
       aes(x = factor(experiment,
                      levels = names(experiment_size)), 
           y = lfc)) + geom_violin() +
  theme_minimal() +
  theme(axis.title = element_text(size=12),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        axis.text = element_text(size=10),
        panel.background = element_rect(fill=NULL,
                                        color = "black"),
        panel.grid = element_blank()) +
  xlab("dataset")

pdf("./analyses/figures/sup/deg_stats.pdf",
    width = 8,
    height = 10)


print(plot_grid(pval,tval,lfc,
          align = "v", ncol = 1,
          rel_heights = c(1,1,1.3)))


dev.off()


