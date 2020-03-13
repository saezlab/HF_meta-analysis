# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de
#' In this script we plot gene level variability,
#' from result objects coming from independent
#' 2-way ANOVAs
#' 

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(tidyverse)
library(cowplot)

anova_confounding = readRDS("data/figure_objects/confounding_anova.rds")

fisher_rank = run_fisher_meta(meta_list = METAheart,
                              n_missing = 4)

genes = names(fisher_rank)

ngenes = 450

# Data frame with ANOVA results

anova_confounding_mod = enframe(lapply(anova_confounding, 
                                       function(x){
                                         colnames(x)[5] = "cfactor"
                                         return(x)
                                       }), "cfactor_name", "values") %>% unnest() %>% 
  dplyr::filter(data_type == "gene_centered") %>% 
  filter(stats=="etasq") %>% 
  mutate(significant = ifelse(gene %in% genes[1:ngenes],
                              "yes","no")) %>% 
  arrange(significant)

# Upper panel

anova_plt_v2 = ggplot(anova_confounding_mod,
                      aes(x=HeartFailure, 
                          y = cfactor, 
                          color = factor(significant,levels = c("no","yes")))) + 
  geom_point(alpha = 2/3) + 
  facet_grid(.~factor(cfactor_name,
                      levels = c("study",
                                 "tech",
                                 "age",
                                 "gender",
                                 "HTx"))) +
  theme_minimal() + scale_color_manual(values=c("lightgrey", "black")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1)) + xlab("eta squared Heart Failure") +
  ylab("eta squared factor")

# Lower panel

anova_boxplt_v2 = ggplot(anova_confounding_mod,
                         aes(y = HeartFailure, 
                             x = factor(significant,levels = c("yes",
                                                               "no")))) + 
  geom_boxplot() + 
  facet_grid(.~factor(cfactor_name,
                      levels = c("study",
                                 "tech",
                                 "age",
                                 "gender",
                                 "HTx"))) +
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        panel.background = element_rect(fill=NULL, 
                                        colour='black',
                                        size=1)) + 
  theme(strip.text.x = element_blank()) +
  xlab("Genes in top 450") +
  ylab("Prop of variance explained by HF")

# Create figure
pdf("analyses/figures/sup/genevar_complete.pdf",
    width = 10,
    height = 8)

plot(plot_grid(anova_plt_v2, anova_boxplt_v2,
               align = "h", ncol = 1,
               rel_heights = c(.5,1)))

dev.off()