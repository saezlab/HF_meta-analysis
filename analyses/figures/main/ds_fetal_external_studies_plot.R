# MIT License

# Copyright (c) [2019] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: This script plots the resutls from the disease score calculation 
# of fetal and external samples (ds_fetal_external_studies.R)

library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork) 
library(gridExtra)
library(grid)


ds_external_experiments = readRDS("data/figure_objects/ds_external_experiment.rds")
ds_external_AUC = readRDS("data/figure_objects/ds_external_experiment_AUC.rds")
cbPalette = c("005073","#d9534f")

#### external plots
# collect information of the AUC performance for each study
s = data.frame("study" = names(ds_external_AUC),"AUC"= NA) %>% column_to_rownames("study")
for (x in names(ds_external_AUC)){
  s[x,1] = ds_external_AUC[[x]]$AUC_All
}

#create tidy data frame for plotting and include AUC into the label
plot.ex.data= enframe(ds_external_experiments, "study") %>% 
  unnest() %>% 
  left_join(s %>% rownames_to_column("study")) %>%
  mutate(AUC = round(AUC,2)) %>%
  mutate(label = paste0(study,'\n',"AUROC: ",AUC))


# 1. Plotting HF of diverse etiologies
plot.ex.data.et = plot.ex.data %>% 
  filter(study !="GSE3586", study != "GSE76701", study != "GSE52601")

# Boxplots of HF of different etiologies
ds.ex.plot.et = ggplot(plot.ex.data.et,
                       aes(x=label, y=Risk_Score, color=HeartFailure)) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed")+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values=cbPalette) + 
  ggtitle("HF Studies with etiological variation") +
  ylim(c(-2.5,2.1))+
  theme_minimal()+
  labs(y= "Disease Score",
       x= "",
       color = "HF")+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_line(size =0.5),
        axis.title.y= element_text(size =13),
        axis.text = element_text(size= 10),
        legend.position = "none")

print(ds.ex.plot.et)

# 2. Plotting HF of technical variation
plot.ex.data.tech = plot.ex.data %>% 
  filter(study == "GSE3586" | study == "GSE52601")

# Boxplots of HF in technical variation
ds.ex.plot.tech = ggplot(plot.ex.data.tech,
                         aes(x=label, y=Risk_Score, color=HeartFailure)) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed")+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values=cbPalette) + 
  ggtitle("HF Studies with technical variation") +
  ylim(c(-2.5,2.1))+
  theme_minimal()+
  labs(y= "",
       x= "",
       color = "HF")+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size= 10),
        axis.line.y = element_line(size = 0.5),
        legend.text = element_text(size=12))
#axis.title.y= element_text(size =13))

print(ds.ex.plot.tech)

##### combine fetal and external plot

ds.ex.fet.plot = plot_grid(ds.ex.plot.et, ds.ex.plot.tech,
                           nrow = 1,
                           rel_widths = c(1.5,1),
                           labels= "AUTO",
                           align = "l")

ds.ex.fet.plot
#Save Figures as PDF

pdf("data/figures/sup/SupplementalFigure12.pdf",
    width = 9,
    height = 4)

plot(ds.ex.fet.plot)


dev.off()


