# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we create the figures of
#' gene centered samples

library(tidyverse)
library(cowplot)
library(gridExtra)
library(grid)

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object
pca_meta_sum_scale = readRDS(file = "data/figure_objects/gcentered_PCA_sum.rds")
pca_plot_df_scale = pca_meta_sum_scale$plot_df

pcs_data = readRDS(file = "data/figure_objects/gcentered_PCs_sum.rds")

# Label order
experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

cbPalette = c("#5e3c58", "#d96459")

pca_plot = ggplot(pca_plot_df_scale, aes(x = PC1, y=PC2, 
                                         color = factor(ExpID,
                                                        levels = names(experiment_size)),
                                         shape = Tech)) + 
  geom_point() + theme_minimal() +
  xlab(paste("PC1",
             as.character(pca_meta_sum_scale$importance[2,1] *100),
             "%")) +
  ylab(paste("PC2",
             as.character(pca_meta_sum_scale$importance[2,2] *100),
             "%")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 1, 0, 0), 'cm'),
        legend.position = c(.1, .4)
        #panel.background = element_rect(fill=NULL, colour='black',
        #                              size=1)
        #axis.line.x = element_line(color="black", size = 1),
        #axis.line.y = element_line(color="black", size = 1)
  ) + labs(color = "Experiment")

pca_legend = get_legend(pca_plot)
pca_plot = pca_plot + theme(legend.position = "none",
                            #axis.line.x = element_line(color="black", size = 0.3),
                            #axis.line.y = element_line(color="black", size = 0.3)
                            panel.background = element_rect(fill=NULL, colour='black',
                                                            size=.3))

box_plot_x = ggplot(pca_plot_df_scale,
                    aes(x = HeartFailure, y = PC1,
                        color = HeartFailure)) + geom_boxplot() +
  theme_minimal() + coord_flip() +
  #c(up,right,bottom,left)
  theme(plot.margin = unit(c(0, 1, 1, 0), 'cm'),
        #   legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        legend.position = c(.01, .8)) +
  ylab(paste("PC1",
             as.character(pca_meta_sum_scale$importance[2,1] *100),
             "%")) +
  scale_colour_manual(values= cbPalette)

hf_legend = get_legend(box_plot_x)
box_plot_x = box_plot_x + theme(legend.position = "none")

box_plot_y = ggplot(pca_plot_df_scale,
                    aes(x = HeartFailure, y = PC2,
                        color = HeartFailure)) + geom_boxplot() +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 0, 0, 1), 'cm'),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  ylab(paste("PC2",
             as.character(pca_meta_sum_scale$importance[2,2] *100),
             "%")) +
  scale_colour_manual(values= cbPalette)

left_panel = plot_grid(box_plot_y,NULL,align = "v",
                       ncol = 1,axis = "b", 
                       rel_heights = c(0.8,0.2))

right_panel = plot_grid(pca_plot,box_plot_x,align = "v",
                        ncol = 1,axis = "b", 
                        rel_heights = c(0.8,0.15))

z_trns_plot = plot_grid(left_panel,right_panel, align = "h",rel_widths = c(0.2,.8))

#title = ggdraw() + draw_label("First two principal components of gene centered combined data sets", fontface='bold')

#plot_grid(title, z_trns_plot, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins

legends = plot_grid(pca_legend,hf_legend,
                    ncol=1,rel_heights =c(1,1)) # rel_heights values control title margins

legends = plot_grid(NULL, legends, nrow = 1, rel_widths = c(.2,1))

final_plot = plot_grid(z_trns_plot,legends,
                       nrow=1, rel_widths=c(1, .25)) # rel_heights values control title margins


plot(final_plot)

## Adding tables

mytheme = gridExtra::ttheme_default(
  core = list(fg_params=list(cex = .8)),
  colhead = list(fg_params=list(cex = .8)))

pcat_plot = tableGrob(pcs_data, rows = NULL, 
                      theme = mytheme)

pdf("./analyses/figures/sup/gcentered_figs.pdf",
    width = 15,
    height = 10)

grid.arrange(final_plot, pcat_plot, ncol = 2, widths = c(3.5,1))

dev.off()





