# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script plots supplementary figure of
# genes_best_performance.R

library(ggplot2)
library(cowplot)

InPerformance_df = readRDS(file = "data/figure_objects/InPerformance_df.rds")
InPerformance_zoom_df = readRDS(file = "data/figure_objects/InPerformance_df_zoom.rds")
OutPerformance_df = readRDS(file = "data/figure_objects/OutPerformance_df.rds")

#Creating plot 
InPerformance_plt = ggplot(InPerformance_df, 
                           aes(x = as.numeric(nin), y = meanAUC, 
                               group = PredictedExperiment, 
                               color = PredictedExperiment)) +
  geom_line() + theme_minimal() + 
  geom_vline(xintercept = 500) + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1)) +
  xlab("Number of top genes included in the signature")

explegend = get_legend(InPerformance_plt)
InPerformance_plt = InPerformance_plt + theme(legend.position = "none")

#Creating plot 

OutPerformance_plt = ggplot(OutPerformance_df, 
                            aes(x = as.numeric(nout), y = meanAUC, 
                                group = PredictedExperiment, 
                                color = PredictedExperiment)) +
  geom_line() + theme_minimal() + 
  geom_vline(xintercept = 500) + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1),
        legend.position = "none") +
  xlab("Number of top genes excluded in the signature")


pdf("./data/figures/sup/SupplementalFigure10.pdf",
    width = 10,
    height = 4.5)

plot(plot_grid(plot_grid(InPerformance_plt, OutPerformance_plt,
          align = "h",
          nrow = 1, rel_widths = c(1,1)),explegend,
          nrow = 1, rel_widths = c(1,.2)))

dev.off()







