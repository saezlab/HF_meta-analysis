# Author: Ricardo Ramirez, 2019
# Description: Figures that correspond to general_variability.R
# Available in supplementary material
#

library(tidyverse)
library(cowplot)

# 1st panel

pca_meta_sum = readRDS(file = "data/figure_objects/pca_meta_summary.rds")

pca_plot_df = pca_meta_sum$plot_df

pca_plot = ggplot(pca_plot_df, aes(x = PC1, y=PC2, 
                                   color = ExpID, shape = Tech)) + 
  geom_point() + theme_minimal() +
  theme(axis.title = element_text(size =12),
        axis.text= element_text(size =12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                      size=1))+
  xlab(paste("PC1",
             as.character(round(pca_meta_sum$importance[2,1] *100)),
             "%")) +
  ylab(paste("PC2",
             as.character(round(pca_meta_sum$importance[2,2] *100)),
             "%"))

# 2nd panel

pca_meta_sum_z = readRDS(file = "data/figure_objects/pca_meta_summary_z.rds")

pca_plot_df = pca_meta_sum_z$plot_df

pca_plot_z = ggplot(pca_plot_df, aes(x = PC1, y=PC2, 
                                   color = ExpID, shape = Tech)) + 
  geom_point() + theme_minimal() +
  theme(axis.title = element_text(size =12),
        axis.text= element_text(size =12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1))+
  xlab(paste("PC1",
             as.character(round(pca_meta_sum_z$importance[2,1] *100)),
             "%")) +
  ylab(paste("PC2",
             as.character(round(pca_meta_sum_z$importance[2,2] *100)),
             "%"))

# 3rd panel

tsne_plotdf = readRDS(file = "data/figure_objects/tsne_z.rds")

tsne_plot_z = ggplot(tsne_plotdf %>% 
                       dplyr::mutate(Study = ExpID),
                     aes(x = tSNE1, y=tSNE2, 
                                     color = Study, shape = Tech)) + 
  geom_point() + theme_minimal() +
  theme(axis.title = element_text(size =12),
        axis.text= element_text(size =12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1))+
  xlab("Dimension 1") +
  ylab("Dimension 2")

# Complete figure

all_legend = get_legend(tsne_plot_z)

pca_plot = pca_plot + theme(legend.position = "none")
pca_plot_z = pca_plot_z + theme(legend.position = "none")
tsne_plot_z = tsne_plot_z + theme(legend.position = "none")

all_panels =plot_grid(pca_plot, pca_plot_z, tsne_plot_z,
            align = "h",
            nrow = 2,
            rel_widths = c(1,1,1),labels = c("A","B","C")
            )

final_plot = plot_grid(all_panels, all_legend, nrow = 1,
          rel_widths = c(1,.2))

# TO DO: define size to export

pdf("./data/figures/sup/SupplementalFigure6.pdf",
    width = 10,
    height = 8)
plot(final_plot)
dev.off()



