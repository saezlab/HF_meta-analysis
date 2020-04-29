# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the lack of difference between DCM and ICM
#' described in the first section of the paper

library(ggplot2)
library(gridExtra)
library(grid)

pca_data = readRDS(file = "data/figure_objects/dcm_icm_pca.rds")
pcs_data = readRDS(file = "data/figure_objects/dcm_icm_pcs.rds")

plotdf = pca_data$plotdf

# 1. Scatterplot

pca_plot = ggplot(plotdf,aes(x = PC1, y=PC2, 
                            color = ExpID, 
                            shape = DCM)) + 
  geom_point() + theme_minimal() +
  theme(axis.title = element_text(size =12),
        axis.text= element_text(size =12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1))+
  xlab(paste("PC1",
             as.character(round(pca_data$importance[2,1] *100)),
             "%")) +
  ylab(paste("PC2",
             as.character(round(pca_data$importance[2,2] *100)),
             "%")) + labs(color = "Study")

# 2. Datatable

mytheme = gridExtra::ttheme_default(
  core = list(fg_params=list(cex = .8)),
  colhead = list(fg_params=list(cex = .8)))

pcat_plot = tableGrob(pcs_data, rows = NULL, 
                      theme = mytheme)


pdf("./data/figures/sup/SupplementalFigure8.pdf",
    width = 13,
    height = 9)

grid.arrange(pca_plot, pcat_plot, ncol = 2, widths = c(3.5,1))


dev.off()
