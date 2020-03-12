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

fisher_rank = run_fisher_meta(meta_list = METAheart ,
                              n_missing = length(METAheart) - 10)

genes = names(fisher_rank)

ngenes = 500

# Scatter plots

anova_confounding_mod = enframe(lapply(anova_confounding, function(x){
  colnames(x)[5] = "cfactor"
  return(x)
}), "cfactor_name", "values") %>% unnest() %>% 
  dplyr::filter(data_type == "gene_centered")


anova_plt_v2 = anova_confounding_mod %>% filter(stats=="etasq") %>% 
  mutate(significant = ifelse(gene %in% genes[1:ngenes],
                              "yes","no")) %>% arrange(significant) %>%
  ggplot(aes(x=HeartFailure, 
             y = cfactor, 
             color = factor(significant,levels = c("no","yes")))) + geom_point(alpha = 2/3) + 
  facet_grid(.~factor(cfactor_name,
                      levels = c("study",
                                 "tech",
                                 "age",
                                 "gender",
                                 "HTx"))) +
  theme_minimal() + scale_color_manual(values=c("lightgrey", "black")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        panel.background = element_rect(fill=NULL, colour='black',
                                        size=1)) + xlab("eta squared Heart Failure") +
  ylab("eta squared factor")


pdf("./analyses/figures/sup/gene_variability_anova.pdf",
    width = 11,
    height = 3)

plot(anova_plt_v2)

dev.off()


















