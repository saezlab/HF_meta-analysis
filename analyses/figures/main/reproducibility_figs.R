# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we generate the replicability 
#' main figure from the manuscript.
#' Tile plots are unified from result objects

library(tidyverse)
library(cowplot)

METAheart = readRDS(file = "data/METAheart.rds") #main object
experiments = names(METAheart)
names(experiments) = experiments

# For labeling
experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

#1. Jaccard plot

jaccard_df = readRDS(file = "data/figure_objects/jaccard_df.rds")


jaccard_tile = ggplot(jaccard_df, aes(x = Var1, 
                                      y = Var2,
                                      fill = value)) +
  geom_tile() +
  theme_minimal() + xlab("Predicted") +
  ylab("Predictor") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        axis.text = element_text(size=10),
        legend.key.size = unit(.6, "cm"),
        legend.text = element_text(size = 9)) +
  scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, 'RdPu'), limits=c(0, .35)) +
  labs(fill = "Jaccard \nIndex")

#2. AUC of Disease scores

pairwise_200 = readRDS(file = "data/figure_objects/pairwise_200.rds")

pcolors = RColorBrewer::brewer.pal(9, 'Purples')[1:6]

pairwise_plot_200 = pairwise_200 %>% ggplot(aes(x=factor(PredictedExperiment,
                                                         levels = rev(names(experiment_size))),
                                                y=factor(PredictorExperiment,
                                                         levels = names(experiment_size))
                                                ,fill = single)) + geom_tile() +
  theme_minimal() + xlab("Predicted") +
  ylab("Predictor") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 9)) +
  scale_fill_gradientn(colours= pcolors, limits=c(0, 1)) + 
  coord_flip() + labs(fill = "AUC")


#3. ES plots

pcolors_up = RColorBrewer::brewer.pal(9, 'Blues')[1:7]
pcolors_down = rev(RColorBrewer::brewer.pal(9, 'Reds')[1:7])
pcolors = c(pcolors_down,pcolors_up)

#upregulation
up_ES = readRDS(file = "data/figure_objects/up_ES.rds")

up_ES_plot = up_ES %>% ggplot(aes(x=factor(Reference,
                                           levels = rev(names(experiment_size))),
                                  y=factor(DEG,
                                           levels = names(experiment_size))
                                  ,fill = ES)) + geom_tile() +
  theme_minimal() + xlab("Reference") +
  ylab("Individual DEG") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9)) +
  scale_fill_gradientn(colours= pcolors, limits=c(-1, 1)) + 
  coord_flip() + ggtitle("Upregulated genes")

up_ES_plot = up_ES_plot + theme(legend.position = "none")

#downregulation

down_ES = readRDS(file = "data/figure_objects/down_ES.rds")

down_ES_plot = down_ES %>% ggplot(aes(x=factor(Reference,
                                           levels = rev(names(experiment_size))),
                                  y=factor(DEG,
                                           levels = names(experiment_size))
                                  ,fill = ES)) + geom_tile() +
  theme_minimal() + xlab("Reference") +
  ylab("Individual DEG") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 9)) +
  scale_fill_gradientn(colours= pcolors, limits=c(-1,1)) + 
  coord_flip() + ggtitle("Downregulated genes")


es_legend = get_legend(down_ES_plot)

down_ES_plot = down_ES_plot + theme(legend.position = "none")


# Align all plots

left_panel = plot_grid(jaccard_tile, pairwise_plot_200,
                       align = "v",
                       ncol = 1, rel_heights = c(1,1),
                       labels = c("A","B"))
right_panel = plot_grid(up_ES_plot,down_ES_plot,
                        align = "v",ncol = 1, 
                        rel_heights = c(1,1),
                        labels = c("C",""))

right_panel = plot_grid(right_panel, es_legend,
                        nrow = 1, rel_widths = c(1,.25))

pdf("./data/figures/main/Figure2.pdf",
    width = 10,
    height = 8)

plot(plot_grid(left_panel, right_panel,
          nrow = 1, rel_widths = c(1,1),
          align = "h"))

dev.off()


























