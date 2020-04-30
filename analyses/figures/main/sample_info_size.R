# MIT License

# Copyright (c) [2019] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: Plotting of Figure 1 (Samplesize, Clinical Characteristics and Gene coverage)

library(tidyverse)
library(ggpubr)
library(ggforce)
library(cowplot)
library(patchwork) 
library(gridExtra)
library(grid)

#loading prerequisites
load("data/clinical_description/dictionaryIDs.ro")
sample = readRDS("data/clinical_description/sample_sizes.rds")%>% 
  as_tibble() %>% 
  mutate(study = as.character(study))
Experiments = dictionary$newID
META  = readRDS("data/METAheart.rds")

clinicalboolean = as.tibble(read.table("data/clinical_description/Tables for MetaHeart Manuscript - Table 2 - Clinical Characteristics.csv", 
                                       header = T, 
                                       sep = ",", 
                                       stringsAsFactors = F)) 

#rename the studies based on the IDs from dictionary
clinicalboolean = clinicalboolean %>% 
  inner_join(.,dictionary %>% rename(Study.ID = GEO_ID)) %>% 
  select(-Study.ID) %>% 
  rename(Study.ID= newID) %>% 
  mutate(Study.ID = factor(Study.ID, levels = c(rev(sample$study))))

##################### Figure1, plotting of table, tile, and sample size

##################### make tile plot for clinical boolean

#Data Transformation
clinic.info = clinicalboolean %>% 
  gather(Age, Sex, Medication, Ethnicity, Comorbidities, EF, NYHA,
         key= "Clinical_Info", 
         value = "Availability") %>% 
  mutate(Clinical_Info = factor(Clinical_Info, levels= c("Age", "Sex", "Medication", "Ethnicity", "Comorbidities", "EF", "NYHA"))) %>%
  mutate(Study.ID = factor(Study.ID, levels = c(rev(sample$study)))) %>%
  mutate(Availability = factor(Availability, levels= c("yes", "no(*)", "no")))

# Data plotting  
plot.tile = ggplot(data= clinic.info, aes(y= Study.ID, x= Clinical_Info, fill = Availability)) +
  geom_tile()+
  scale_fill_manual(values = c("#444c5c","#b0aac0","#f0efef"))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, 
                                   hjust = 0),
        axis.title.x.top = element_text(),
        axis.title.y = element_blank(),
        text= element_text(size= 15))+
  labs(x="")+
  scale_x_discrete(position = "top") 

plot.tile

##################### Sample size 

sample.size = readRDS(file="data/clinical_description/sample_sizes.rds") 

sample.plot = sample.size %>% 
  gather(key = Sample, value = number , CT, ICM, DCM)  %>%
  mutate(study= factor(study,  levels = rev(as.character(sample.size$study))))%>%
  mutate(Sample= factor(Sample, levels = c("ICM", "DCM", "CT")))

# Data plotting
plot.samplesize = ggplot(data = sample.plot, 
                         aes(x= study, 
                             y= number, 
                             fill = Sample))+
  geom_col()+
  scale_fill_manual(values = c("#2b9679","#54ce93","#1e4d62"))+ 
  coord_flip()+
  labs(x= "",
       y= "Sample size",
       fill = "Sampletype")+
  theme_minimal()+
  scale_y_continuous(position = "bottom") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 15),
        axis.title.x= element_text(size= 12))

plot.samplesize
        
##################### Combine Tile and Sample plots ##################### 

leg1 <- get_legend(plot.tile + 
                     guides(color = guide_legend(nrow = 1)) +
                     theme(legend.position = "bottom"))
leg2 <- get_legend(plot.samplesize+ 
                     guides(color = guide_legend(nrow = 1)) +
                     theme(legend.position = "bottom"))

plot.samplesize = plot.samplesize +theme(legend.position = "none")
plot.tile = plot.tile  +theme(legend.position = "none")

plot.tilesample = plot_grid(plot.tile, plot.samplesize,
                      ncol = 2,
                      align = "h",
                      labels = "AUTO", 
                      rel_widths = c(1,1.5))

blank_p <- plot_spacer() + theme_void()

leg12 <- plot_grid(blank_p,leg1,blank_p, leg2,
                   blank_p, blank_p,
                   nrow = 1)

fig_1 = plot_grid(plot.tilesample, leg12, nrow = 2,
                    align = "h",
                    axis = "t",
                    rel_widths = c(1, 0.1), 
                    rel_heights = c(1,0.1))
fig_1

####################### print PDFs

pdf("data/figures/main/Figure1.pdf",
    width = 10,
    height = 6)
fig_1
dev.off()
