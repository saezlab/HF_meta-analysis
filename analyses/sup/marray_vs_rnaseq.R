# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare microarrays vs RNASeq
#' 


source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
library(cowplot)

# Load dictionary of studies
METAheart = readRDS(file = "data/METAheart.rds") #main object

# Get marrays IDs
marrays = c("GSE76701","GSE57345","GSE42955",
            "GSE1869","GSE3585","GSE26887",
            "GSE5406","GSE16499")

load("./data/dictionaryIDs.ro")
dictionary = dictionary %>% 
  dplyr::mutate(tech = ifelse(GEO_ID %in% marrays,
                              "marray","rnaseq"))

# Define objects

# Compare gene coverage

gcoverage = map(METAheart, function(x){
  return(nrow(x$GEX))
}) %>% enframe("study",value = "gene_coverage") %>% 
  unnest() %>% left_join(dictionary, by = c("study" = "newID"))

ggplot(gcoverage, aes(x = tech, y = gene_coverage, label = study)) + geom_violin(draw_quantiles = c(0.25,.5,.75)) +
  geom_point(aes(color = study)) + ggrepel::geom_text_repel() + theme_minimal()

# Compare gene coverage in ReHeat

fisher_rank = readRDS(file = "data/shiny/fisher_rank.rds")

gcoverage_meta = map(METAheart, function(x){
  return(length(intersect(rownames(x$GEX),names(fisher_rank))))
}) %>% enframe("study",value = "gene_coverage") %>% 
  unnest() %>% left_join(dictionary, by = c("study" = "newID"))


# Compare sample size
ssize = map(METAheart, function(x){
  return(ncol(x$GEX))
}) %>% enframe("study",value = "sample_size") %>% 
  unnest() %>% left_join(dictionary, by = c("study" = "newID"))

#ggplot(ssize, aes(x = tech, y = sample_size,
#                  label = study)) + geom_violin(draw_quantiles = c(0.25,.5,.75)) +
#  geom_point(aes(color = study)) + ggrepel::geom_text_repel() + theme_minimal()

# Compare variance captured by HF 

HFvariance = readRDS(file = "data/figure_objects/covariate_var_all.rds") %>%
  left_join(dictionary, by =c("study"="newID")) %>%
  dplyr::filter(name == "HeartFailure")
  
#ggplot(HFvariance, aes(x = tech, 
#                      y = total_var,
#                      label = study)) + 
#  geom_violin(draw_quantiles = c(0.25,.5,.75)) +
#  geom_point(aes(color = study)) + ggrepel::geom_text_repel() + theme_minimal()

# "Contribution" to the current ranking

contribution = readRDS("data/figure_objects/contribution.rds") %>%
  dplyr::select(pathway, ES) %>% 
  left_join(dictionary, by = c("pathway" = "newID"))

#ggplot(contribution, aes(x = tech, 
#                       y = ES,
#                       label = pathway)) + 
#  geom_violin(draw_quantiles = c(0.25,.5,.75)) +
#  geom_point(aes(color = pathway)) + ggrepel::geom_text_repel() + theme_minimal()

# Predictiveness of single-study
# label misleading, this is top 500

AUC_performance = readRDS("data/figure_objects/pairwise_200.rds") %>%
  dplyr::filter(PredictorExperiment != PredictedExperiment) %>%
  dplyr::mutate(Experiment = PredictorExperiment) %>%
  group_by(Experiment) %>%
  dplyr::summarise("Predictor_performance" = mean(single)) %>%
  left_join(readRDS("data/figure_objects/pairwise_200.rds") %>%
              dplyr::filter(PredictorExperiment != PredictedExperiment) %>%
              dplyr::mutate(Experiment = PredictedExperiment) %>%
              group_by(Experiment) %>%
              dplyr::summarise("Predicted_performance" = mean(single)))

# Jaccard Index
# label misleading, this is top 500

jaccard_ix = readRDS(file = "data/figure_objects/jaccard_df.rds") %>%
  dplyr::filter(Var1 != Var2)

# lou results

jaccard_LOU_rankings = readRDS(file = "data/figure_objects/jaccard_LOU_rankings.rds") %>%
  dplyr::filter(StudyA != StudyB) %>%
  group_by(StudyA) %>%
  summarise(mean_Jaccard_lou = mean(JaccardIx))

# All data together

all_chars = left_join(left_join(left_join(left_join(left_join(gcoverage, ssize),HFvariance), contribution),
                      AUC_performance, by = c("study"="Experiment")),
                      jaccard_LOU_rankings, by = c("study"="StudyA"))

# Supplemental plots

#Wilcoxon tests 
wilcox.test(ES ~ tech, data = all_chars)
wilcox.test(Predictor_performance ~ tech, data = all_chars)

upper_panel = all_chars %>%
  dplyr::select(study, tech, ES, Predictor_performance) %>%
  pivot_longer(c(ES,Predictor_performance),names_to = "name",
               values_to = "value") %>%
  mutate(name = ifelse(name == "ES",
                        "ES of individual DEGs in the top of HF-CS",
                        "mean AUROC of prediction of HF labels in rest")) %>%
  mutate(name = factor(name,
                       levels =c("mean AUROC of prediction of HF labels in rest",
                       "ES of individual DEGs in the top of HF-CS"))) %>%
  ggplot(aes(x = tech, y = value,label = study)) + 
  geom_boxplot() + geom_point(aes(color = study)) + 
   ggrepel::geom_text_repel(size = 3) + theme_classic() +
  facet_grid(.~name,scales = "free") +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill=NULL,
                                        color = "black"),
        panel.grid = element_blank())

# correlations

# No correlation between HF variability and predictability
cor.test(all_chars$total_var,
         all_chars$Predictor_performance,method = "spearman")

cor.test(all_chars$sample_size,
         all_chars$Predictor_performance,method = "spearman")

lower_panel = all_chars %>%
  dplyr::select(study, tech, Predictor_performance,
                sample_size,total_var) %>%
  pivot_longer(c(sample_size,total_var),names_to = "name",
               values_to = "value") %>%
  mutate(name = ifelse(name == "sample_size",
                       "Sample Size",
                       "Proportion of variance explained by HF")) %>%
  ggplot(aes(x = value, y = Predictor_performance,label = study)) + 
  geom_point(aes(color = study)) + 
  ggrepel::geom_text_repel(size = 3) + theme_classic() +
  facet_grid(.~name,scales = "free") +
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill=NULL,
                                        color = "black"),
        panel.grid = element_blank()) +
  ylab("Prediction performance (mean AUROC)")

pdf(file = "data/figures/sup/marrayvsrnaseq.pdf", 
    width = 8, height = 7)

plot_grid(upper_panel,lower_panel,ncol = 1, labels = c("A","B"))

dev.off()
