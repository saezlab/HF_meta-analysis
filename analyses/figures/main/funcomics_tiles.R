# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the results of funcomics 
#' applied to the meta-ranked gene list.
#' Main figure
#' 
#' 

library(tidyverse)
library(cowplot)

GSEA = readRDS(file = "data/figure_objects/GSEA_results.rds")
dorothea = readRDS(file = "data/figure_objects/dorothea_results.rds")
progeny = readRDS(file = "data/figure_objects/PROGENy_results.rds")
miRNA = readRDS(file = "data/figure_objects/GSEA_mir_results.rds")

# aesthetics
get_lfc_colors = circlize::colorRamp2(seq(-4,4,.5),
                                      c(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(9)),
                                        colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds'))(8)))

cols = get_lfc_colors(seq(-4,4,1))

# GSEA filtering

GSEA_filtering = GSEA %>% dplyr::filter(pval < 0.05) %>%
  slice(1:50) %>%
  arrange(-log2(pval) * sign(NES)) %>%
  mutate(dir_reg = sign(ES),
         log_pval = -log10(pval),
         dir_col = dir_reg * log_pval,
         gset = pathway,
         tile_color = get_lfc_colors(dir_col)) %>%
  dplyr::select(gset,dir_col,tile_color)

GSEA_filtering$gset = factor(GSEA_filtering$gset,
                             levels = GSEA_filtering$gset)

g1 <- GSEA_filtering %>% 
  ggplot(aes(x=gset, y=1, fill=factor(1:nrow(GSEA_filtering)),
             height=.5)) + 
  geom_tile() + theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=as.character(GSEA_filtering$tile_color), 
                    na.value='black', guide=FALSE) +
  coord_flip()

# Dorothea_filtering

dorothea_filtering = dorothea %>% dplyr::filter(adj_pvalue < 0.1) %>%
  arrange(-log2(pvalue) * sign(NES)) %>%
  mutate(dir_reg = sign(NES),
         log_pval = -log10(pvalue),
         dir_col = dir_reg * log_pval,
         gset = RegulonName,
         tile_color = get_lfc_colors(dir_col))

dorothea_filtering$gset = factor(dorothea_filtering$gset,
                                 levels = dorothea_filtering$gset)

g2 <- dorothea_filtering %>% 
  ggplot(aes(x=gset, y=1, 
             fill=factor(1:nrow(dorothea_filtering)),
             height=.5)) + 
  geom_tile() + theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 8.5, 0.5), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=as.character(dorothea_filtering$tile_color), 
                    na.value='black', guide=FALSE) +
  coord_flip()

plot_grid(g1,g2)


# PROGENy filtering

progeny_filtering = progeny %>%
  arrange(-log2(progeny_pvals) * sign(progeny_scores)) %>%
  mutate(dir_reg = sign(progeny_scores),
         log_pval = -log10(progeny_pvals),
         dir_col = dir_reg * log_pval,
         gset = pathway,
         tile_color = get_lfc_colors(dir_col))

progeny_filtering$gset = factor(progeny_filtering$gset,
                                 levels = progeny_filtering$gset)

g3 = progeny_filtering %>% 
  ggplot(aes(x=gset, y=1, 
             fill=factor(1:nrow(progeny_filtering)),
             height=.5)) + 
  geom_tile() + theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 12.7, .5), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=as.character(progeny_filtering$tile_color), 
                    na.value='black', guide=FALSE) +
  coord_flip()

# miRNAs

miRNA_filtering = miRNA %>% dplyr::filter(pval < 0.05) %>%
  arrange(-log2(pval) * sign(NES)) %>%
  mutate(dir_reg = sign(ES),
         log_pval = -log10(pval),
         dir_col = dir_reg * log_pval,
         gset = pathway,
         tile_color = get_lfc_colors(dir_col)) %>%
  dplyr::select(gset,dir_col,tile_color)

miRNA_filtering$gset = unlist(lapply(strsplit(miRNA_filtering$gset,
                                              "[A-Z]_"), 
                              function(x) x[2]))

miRNA_filtering$gset = factor(miRNA_filtering$gset,
                             levels = miRNA_filtering$gset)

g4 <- miRNA_filtering %>% 
  ggplot(aes(x=gset, y=1, fill=factor(1:nrow(miRNA_filtering)),
             height=.5)) + 
  geom_tile() + theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 13, .5), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=as.character(miRNA_filtering$tile_color), 
                    na.value='black', guide=FALSE) +
  coord_flip()

pdf("./analyses/figures/main/funcomics_tiles.pdf",
    width = 7.51,
    height = 7.21)

plot(plot_grid(g1,g2,g3,g4,nrow = 1,
          rel_widths = c(1.05,.35,.35,.4)))

dev.off()


## Make legend

legend_df = tibble(size_lab = seq(4,-4,-.5),
                   tile_color = get_lfc_colors(size_lab),
                   logpval = c(seq(4,.5,-.5),0,
                               seq(.5,4,.5))) %>%
  arrange(size_lab)

legend_df$size_lab = factor(legend_df$size_lab,
                            levels = legend_df$size_lab)

legend_plot = legend_df %>% 
  ggplot(aes(x=size_lab, y=1, fill=factor(1:nrow(legend_df)),
             height=.5)) + 
  geom_tile() + theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(13, 10, 0, .5), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_discrete(breaks = seq(4,-4,-2),
                   labels = as.character(c(seq(4,0,-2),
                                           seq(2,4,2)))) + 
  scale_fill_manual(values=as.character(legend_df$tile_color), 
                    na.value='black', guide=FALSE) +
  coord_flip()

pdf("./analyses/figures/main/funcomics_tiles_legend.pdf",
    width = 4.76,
    height = 7.21)

plot(legend_plot)

dev.off()

