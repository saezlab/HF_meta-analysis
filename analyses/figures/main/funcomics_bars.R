# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the results of funcomics 
#' applied to the meta-ranked gene list.
#' Main figure
#' 
#' 

library(tidyverse)
library(cowplot)

GSEA = readRDS(file = "data/shiny/GSEA_results.rds")
dorothea = readRDS(file = "data/figure_objects/dorothea_results.rds")
progeny = readRDS(file = "data/figure_objects/PROGENy_results.rds")
miRNA = readRDS(file = "data/figure_objects/GSEA_mir_results.rds")

# GSEA filtering

GSEA_filtering = GSEA %>% 
  arrange(desc(abs(NES))) %>%
  slice(1:50) %>%
  arrange(-log10(pval) * sign(NES)) %>%
  mutate(pathway = strtrim(pathway,40)) %>%
  mutate(dir_reg = sign(ES),
         log_pval = -log10(pval),
         gset = pathway) %>%
  mutate(gset = ifelse(padj<=0.15,
                       paste("*",gset,sep=""),
                       gset))

GSEA_filtering$gset = factor(GSEA_filtering$gset,
                             levels = GSEA_filtering$gset)

g1 = ggplot(GSEA_filtering, aes(x=gset, y=-log10(pval),
                                 fill = factor(dir_reg))) + 
  geom_bar(stat = "identity") + theme_minimal() + 
  geom_hline(yintercept = -log10(0.05),linetype="dotted" ) +
  coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 9),
        #plot.margin = unit(c(0, 0, 0, 0), 'cm')
        ) + ylab("-log10(p-value)") +
  scale_fill_manual(values = c("steelblue","indianred"))
  

# Dorothea_filtering

dorothea_filtering = dorothea %>% dplyr::slice(1:50) %>%
  #dplyr::filter(pvalue < 0.05) %>%
  arrange(-log10(pvalue) * sign(NES)) %>%
  mutate(dir_reg = sign(NES),
         log_pval = -log10(pvalue),
         dir_col = dir_reg * log_pval,
         gset = RegulonName) %>%
    mutate(gset = ifelse(adj_pvalue<=0.15,
                         paste("*",gset,sep=""),
                         gset))

dorothea_filtering$gset = factor(dorothea_filtering$gset,
                                 levels = dorothea_filtering$gset)

g2 <- ggplot(dorothea_filtering, aes(x=gset, y=-log10(pvalue),
                                 fill = factor(sign(NES)))) + 
  geom_bar(stat = "identity") + theme_minimal() +
  geom_hline(yintercept = -log10(0.05),linetype="dotted" ) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 1),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 9),
        legend.position = "none"
        #plot.margin = unit(c(0, 0, 0, 0), 'cm')
        ) +
  scale_fill_manual(values = c("steelblue","indianred")) +
  coord_flip() + ylab("-log10(p-value)")

# miRNAs

miRNA_filtering = miRNA  %>% dplyr::slice(1:50) %>%
  #dplyr::filter(pvalue < 0.05) %>%
  arrange(-log10(pvalue) * sign(NES)) %>%
  mutate(dir_reg = sign(NES),
         log_pval = -log10(pvalue),
         dir_col = dir_reg * log_pval,
         gset = RegulonName)

miRNA_filtering$gset = unlist(lapply(strsplit(miRNA_filtering$gset,
                                              "[A-Z]_"), 
                                     function(x) x[2]))

miRNA_filtering$gset = unlist(lapply(strsplit(miRNA_filtering$gset,
                                              "_"), 
                                     function(x) x[1]))

miRNA_filtering = mutate(miRNA_filtering, gset = ifelse(adj_pvalue<=0.15,
                     paste("*",gset,sep=""),
                     gset))

miRNA_filtering$gset = factor(miRNA_filtering$gset,
                              levels = miRNA_filtering$gset)

g3 <- ggplot(miRNA_filtering, aes(x=gset, y=-log10(pvalue),
                                     fill = factor(sign(NES)))) + 
  geom_bar(stat = "identity") + theme_minimal() +
  geom_hline(yintercept = -log10(0.05),linetype="dotted" ) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.text.y = element_text(hjust = 1),
        panel.grid = element_blank(),
        legend.position = "none"
        #plot.margin = unit(c(0, 0, 0, 0), 'cm')
  ) +
  scale_fill_manual(values = c("steelblue","indianred")) +
  coord_flip() + ylab("-log10(p-value)")

upper_panel = plot_grid(g1,g2,g3,align = "h", rel_widths = c(1,.45,.45),
          nrow = 1)

# PROGENy filtering

progeny_filtering = progeny %>%
  arrange(-log10(progeny_pvals + 0.00000001) * sign(progeny_scores)) %>%
  mutate(dir_reg = sign(progeny_scores),
         log_pval = -log10(progeny_pvals + 0.00000001),
         gset = pathway,
         adj_pval = p.adjust(progeny_pvals,method = "BH")) %>%
  mutate(gset = ifelse(adj_pval<=0.15,
                       paste("*",gset,sep=""),
                       gset)) %>%
  arrange((progeny_scores))

progeny_filtering$gset = factor(progeny_filtering$gset,
                                levels = progeny_filtering$gset)

g4 = ggplot(progeny_filtering,
            aes(x=gset, y = log_pval,
                fill = factor(dir_reg))) + 
  geom_bar(stat = "identity") + theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_text(size=10),
        legend.direction='vertical',
        plot.margin = unit(c(0, 0, 8, 0), 'cm')
  ) +
  scale_fill_manual(values = c("steelblue","indianred")) + 
  ylab("-log10(p-value)") + 
  geom_hline(yintercept = -log10(0.05),linetype="dotted" ) +
  coord_flip() + labs(fill="Direction")


upper_panel = plot_grid(g1,g2,g3,g4,axis = "t" ,align = "h", rel_widths = c(.7,.35,.30,.30),
                        nrow = 1)
# Final panel

pdf("./data/figures/main/Figure4_raw.pdf",
    width = 10.3,
    height = 7.21)

plot(upper_panel)

dev.off()
