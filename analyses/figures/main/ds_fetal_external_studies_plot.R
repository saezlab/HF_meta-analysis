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

# load output from ds_fetal_external_studies.R)
ds_external_experiments = readRDS("data/ds_external_experiment.rds")
ds_external_AUC = readRDS("data/ds_external_experiment_AUC.rds")
ds_fetal_experiments = readRDS("data/ds_fetal_experiment.rds")
ds_fetal_AUC = readRDS("data/ds_fetal_experiment_AUC.rds")

# colour definition
cbPalette = c("005073","#d9534f")


#### Plotting Figure 5.
# Title: Disease score calculation based on the top 500 genes from the consensus signature 
# for diverse HF studies.

# collect information of the AUC performance for each study
auc = data.frame("study" = names(ds_external_AUC),"AUC"= NA) %>% column_to_rownames("study")
for (x in names(ds_external_AUC)){
auc[x,1] = ds_external_AUC[[x]]$AUC_All
  }

#create tidy data frame for plotting and include AUC into the label
plot.ex.data= enframe(ds_external_experiments, "study") %>% 
  unnest() %>% 
  left_join(auc %>% rownames_to_column("study")) %>%
  mutate(AUC = round(AUC,2)) %>%
  mutate(label = paste0(study,'\n',"AUROC: ",AUC))

# filter for studies with etiological variation
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
  ggtitle("HF studies with diverse etiologies") +
  ylim(c(-2.5,2.1))+
  theme_minimal()+
  labs(y= "Disease score",
       x= "",
       color = "HF")+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_line(size =0.5),
        axis.title.y= element_text(size =13),
        axis.text = element_text(size= 11.5),
        legend.position = "none")

  print(ds.ex.plot.et)

  # filter HF studies with technical variation
  plot.ex.data.tech = plot.ex.data %>% 
    filter(study == "GSE3586" | study == "GSE52601" | study == "GSE76701")
  
  # Boxplots of HF in technical variation
  ds.ex.plot.tech = ggplot(plot.ex.data.tech,
                         aes(x=label, y=Risk_Score, color=HeartFailure)) +
    geom_hline(yintercept = 0,
               color = "grey",
               linetype = "dashed")+
    geom_boxplot(position=position_dodge(0.8))+
    geom_jitter(position=position_dodge(0.8)) +
    scale_colour_manual(values=cbPalette, 
                        name = "",
                        labels = c("CT", "HF")) + 
    ggtitle("HF studies with technical variation") +
    ylim(c(-2.5,2.1))+
    theme_minimal()+
    labs(y= "",
         x= "",
         color = "HF")+
    theme(panel.grid.major = element_blank(),
          axis.text = element_text(size= 11.5),
          axis.line.y = element_line(size = 0.5))
          #axis.title.y= element_text(size =13))
  
  print(ds.ex.plot.tech)
  
  
## combine previous plots for Figure 5
 ds.ex.plot = plot_grid(ds.ex.plot.et, ds.ex.plot.tech, 
                             nrow = 1,
                             rel_widths = c(1.2,1),
                             labels= "AUTO",
                             align = "l")
  
#Save Figures as PDF
  pdf("figures_pdfs/ds_external_studies.pdf",
      width = 10,
      height = 4)
  
  ds.ex.plot
  
  dev.off()
  
    
#### Plotting of the disease score classification of fetal samples
  
  # collect information of the AUC performance for each study
  auc = data.frame("study" = names(ds_fetal_AUC),"AUC"= NA) %>% column_to_rownames("study")
  for (x in names(ds_fetal_AUC)){
    auc[x,1] = ds_fetal_AUC[[x]]$AUC_All
  }
  
  #create tidy data frame for plotting and include AUC into the label
  plot.fet.data= enframe(ds_fetal_experiments, "study") %>% 
    unnest() %>% 
    left_join(auc %>% rownames_to_column("study")) %>%
    mutate(AUC = round(AUC,2)) %>%
    mutate(label = paste0(study,'\n',"AUROC: ",AUC))

#plot.fet.data$label[grepl("PRJNA522417", plot.fet.data$label)] = "Spurrell19"
plot.fet.data$label= gsub("PRJNA522417", "Spurrell19", plot.fet.data$label)

# Boxplots of fetal datasets
ds.fet.plot = ggplot(plot.fet.data,
                     aes(x=label, y=Risk_Score, color=HeartFailure)) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed")+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values=c("005073","#499483"),
                      name = "",
                      labels = c("CT", "fetal") ) + 
  labs(y= "Disease Score",
       x= "",
       color = "Fetal Sample")+
  theme_minimal()+
  ylim(c(-2.5,2.1))+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_line(size =0.5),
        axis.text = element_text(size= 11.5),
        axis.title.y= element_text(size =13))
        
        print(ds.fet.plot)


#Save Figures as PDF

pdf("figures_pdfs/ds_fetal.pdf",
    width = 4,
    height = 4)

ds.fet.plot

dev.off()


