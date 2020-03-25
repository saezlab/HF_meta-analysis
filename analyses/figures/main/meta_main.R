# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: This script plots the main figure of our meta-analysis.
# For the top most consistent genes, we plot their experiment specific
# LFCs. We add the ranking p-value to show how the ranking was built.
# Code used as template comes from:
# [https://github.com/zellerlab/crc_meta]

library(ggplot2)
library(cowplot)

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

marker_genes = c("MYH6","MME","CNN1","NPPA","KCNH2","SLC2A1",
                 "ATP2A2","COL21A1","COL15A1","ECM2","MXRA5",
                 "KIT","FNDC1","LAMA4","SSPN","KCNN3","FGF14")

experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

experiments = experiments[names(experiment_size)]

# aesthetic parameters: Generate LFC color function generator and legend

get_lfc_colors = circlize::colorRamp2(seq(-3,3,0.25),
                                      c(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(12)),
                                        "white",
                                        colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds'))(12))
)

axis_text_lfc = seq(-5,5,1)

lfc_legend_df = tibble(LFC=seq(-6,6,.1)) %>%
  mutate(lfc_color = get_lfc_colors(LFC)) %>%
  mutate(axis_text = (ifelse(LFC %in% axis_text_lfc,
                             as.character(LFC),""))) %>%
  mutate(LFC = factor(LFC, levels = LFC))

lfc_legend_plot = lfc_legend_df %>% 
  ggplot(aes(x=LFC, y=1, fill=factor(1:nrow(lfc_legend_df)))) + 
  geom_tile() + theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 15),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=lfc_legend_df$lfc_color, 
                    na.value='black', guide=FALSE,
                    breaks =  c(seq(3,-3,-1.5))) + 
  scale_x_discrete(labels=lfc_legend_df$axis_text) + 
  ggtitle("LFC") + coord_flip()

pdf("analyses/figures/main/meta_analysis_hmps_simple_legend.pdf",
    width = 1,
    height = 4)

plot(lfc_legend_plot)

dev.off()

# meta analysis
fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = length(METAheart) - 10) #Fisher combined test

# Choosing useful genes
topgenes = names(fisher_rank)[1:500]

t_matrix = get_all_limma(meta_list = METAheart,
                         limma_column = "logFC")

t_rank = sort(rowMeans(t_matrix[topgenes,],na.rm = T),
              decreasing = F)

# Processing meta-heart object
meta_limma_summary = lapply(METAheart , function(x, genes){
  limma_results = x[["HF_limma"]] 
  ix = fastmatch::fmatch(genes, limma_results[["ID"]]) #Get genes ix
  
  limma_results_tibble = tibble(id = factor(genes,levels = genes),
                                t = limma_results[["t"]][ix],
                                logfc = limma_results[["logFC"]][ix],
                                p_val = limma_results[["P.Value"]][ix]) %>%
    mutate(lfc_color = get_lfc_colors(logfc))
  
  return(limma_results_tibble)
}, genes = names(t_rank))

# 3.5. Make cowplot of hmaps

plot_limma_summary = function(plot_tibble,study_name){
  g1 = plot_tibble %>% 
    ggplot(aes(x=id, y=1, fill=factor(1:nrow(plot_tibble)))) + 
    geom_tile() + theme_minimal() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          #axis.title = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 360,vjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='white'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          legend.position = "none") + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_manual(values=as.character(plot_tibble$lfc_color), 
                      na.value='black', guide=FALSE) +
    #ggtitle(study_name)
    ylab(study_name)
  return(g1)
}


plot_list = lapply(experiments,function(x){
  plot_tibble = meta_limma_summary[x]
  plot_limma_summary(plot_tibble[[1]], study_name = x)
})

# Last plot 
#highlight = sample(topgenes,5)
last_tibble = meta_limma_summary[experiments[13]][[1]] %>%
  mutate(lfc_color = ifelse(id %in% marker_genes,
                            "white","white"))
# Here define the genes you want
last_g = last_tibble %>% 
  ggplot(aes(x=id, y=1, fill=factor(1:nrow(last_tibble)))) + 
  geom_tile() + theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 30, 
                                   hjust = 1,
                                   size=7), 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='white'),
        plot.margin = unit(c(0, 0, 0, .5), 'cm'),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_discrete(breaks = marker_genes) + #Here add breaks with genes
  scale_fill_manual(values=as.character(last_tibble$lfc_color), 
                    na.value='black', guide=FALSE)


gene_expression_plot = plot_grid(plot_list[[1]],plot_list[[2]], plot_list[[3]],
                                 plot_list[[4]],plot_list[[5]], plot_list[[6]],
                                 plot_list[[7]],plot_list[[8]], plot_list[[9]],
                                 plot_list[[10]],plot_list[[11]], plot_list[[12]],
                                 plot_list[[13]],plot_list[[14]],plot_list[[15]],
                                 plot_list[[16]],
                                 last_g,
                                 ncol = 1, align = 'v', rel_widths = c(rep(0.12, 13),0.34))

pdf("analyses/figures/main/meta_analysis_hmps_simple.pdf",
    width = 8.70,
    height = 3.54)

plot(gene_expression_plot)

dev.off()

# p_value plot

pvalue_df = tibble(p_val = -log10(fisher_rank),
                   rank = 1:length(p_val))

rank_plot = ggplot(pvalue_df, aes(x = rank,
                      y = p_val)) +
  geom_point() + 
  theme_minimal() +
  theme(axis.text = element_text(size=10)) +
  scale_x_continuous(expand = c(0.01, 0),breaks = c(length(topgenes),
                                                    seq(1000,3000,2000), 
                                               seq(3000,11000,2000),
                                               length(fisher_rank)),
                                    limits = c(0.3,14100)) +
  geom_vline(xintercept = 500) +
  ylab("-log10(BH p-value)") +
  xlab("Meta-analysis Rank")


pdf("analyses/figures/main/meta_analysis_rank.pdf",
    width = 8.70,
    height = 4.24)

plot(rank_plot)

dev.off()







































