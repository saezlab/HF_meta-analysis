# MIT License

# Copyright (c) [2019] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: Plotting of Gender and Age pre HF Study in MetaHeart project

library(tidyverse)
library(ggpubr)
library(ggforce)
library(formattable)
library(cowplot)
library(patchwork) 
library(gridExtra)
library(grid)

#load prerequisites
load("data/clinical_description/dictionaryIDs.ro")
sample = readRDS("data/clinical_description/sample_sizes.rds")%>% 
  as_tibble() %>% 
  mutate(study = as.character(study))
Experiments = dictionary$newID

clinicalboolean = as.tibble(read.table("data/clinical_description/Tables for MetaHeart Manuscript - Table 2 - Clinical Characteristics.csv", 
                                       header = T, 
                                       sep = ",", 
                                       stringsAsFactors = F)) 


#rename the studies based on the IDs from dictionary
clinicalboolean = clinicalboolean %>% 
  inner_join(.,dictionary %>% rename(Study.ID = GEO_ID)) %>% 
  dplyr::select(-Study.ID) %>% 
  rename(Study.ID= newID) %>% 
  mutate(Study.ID = factor(Study.ID, levels = c(rev(sample$study))))

##################### Plot Ages ##################### 
### Data preparation

load("data/clinical_description/ClinicalCharacteristics_age.ro")

# transform ages data frame into a tidy data frame
ages_tidy = ages %>% rownames_to_column("Study") %>% 
  as.tibble()%>% 
  rename("HF" = "HFMeanAge") %>% 
  rename("CT" = "NFMeanAge") %>%
  gather(HF, CT, key = "Group", value= "Age") %>%
  rename("HF" = "HFSDAge") %>% 
  rename("CT"= "NFSDAge") 

agesd = ages_tidy %>% 
  dplyr::select(Study, HF, CT)

agesd = agesd[1:16,] %>%
  gather(HF, CT, value = AgeSD, key= Group)

ages_tidy = agesd %>% 
  inner_join(ages_tidy %>% 
               select(Group, Age, Study), by = c("Group", "Study"))

# alphamatrix contains summarystatistic information
alphamatrix = clinicalboolean %>% 
  select(Study.ID, Age) %>% 
  rename(Study = Study.ID) %>%
  rename(Samplestatistic=Age)

# agesplot is the final df for plotting, containing also the summarystatistic information
agesplot = ages_tidy %>% 
  left_join(alphamatrix, by= "Study") %>%
  mutate(Study = factor(Study, levels = c(sample$study)))

# summarystatistic information gets formatted, to be used for the alpha level in plots 
agesplot$Samplestatistic= gsub("[(*)]", "" , agesplot$Samplestatistic)
agesplot$Samplestatistic= gsub("yes", "available", agesplot$Samplestatistic)
agesplot$Samplestatistic =gsub("no", "not available",agesplot$Samplestatistic)


### Data plotting
plot.Age = ggplot(agesplot)+
  facet_grid(. ~ Study)+
  geom_col(aes(y=Age,
               x=Group, 
               #alpha = Samplestatistic,
               fill = Group))+
  scale_alpha_discrete(range = c(1,0.5))+
  #scale_fill_manual(values = c("#444c5c","#499483"))+
  scale_fill_manual(values =c("#7c5074","#7c5074"))+
  geom_errorbar(aes(x = Group, 
                    ymin=Age-AgeSD, 
                    ymax=Age+AgeSD), 
                width=.2)+
  theme_minimal()+
  labs(x="", 
       y = "Age (y)", 
       alpha = "Sample Information",
       fill= "Disease")+
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #axis.text.x=element_blank(),
        legend.position = "none",
        strip.text = element_text(angle = 90,
                                  vjust = 0,
                                  hjust= 0),
        axis.text.x=element_blank(),
        text = element_text(size=15))+
  coord_cartesian(ylim=c(22,73))

plot.Age

##################### Plot Gender ##################### 
### Data preparation
load("data/clinical_description/ClinicalCharacteristics_gender.ro")

#transform gender data frame into a tidy data frame
gendertidy = gender %>% 
  rownames_to_column("Study") %>% 
  as.tibble() %>% 
  rename(HF = female_HF) %>% 
  rename(CT= female_NF) %>% 
  gather(HF, CT, value= Female, key = Group) %>%
  mutate(Male = 100-Female) %>% 
  gather(Female, Male, value = Percentage, key= Gender)

#aplhamatrix is the information for the alpha value (if a study provides samplewise information or only summary statistics)
alphamatrix = clinicalboolean %>% 
  as.tibble() %>%
  select(Study.ID, Sex) %>% 
  rename(Study= Study.ID) %>%
  rename(Samplestatistic = Sex)

# genderplot is the joined alphamatrix and gendertidy, i.e. the data frame in tidy format to plot
genderplot = gendertidy %>% 
  left_join(alphamatrix) %>%
  mutate(Study = factor(Study, levels = c(sample$study)))


### Data plotting
plot.Gender = ggplot(genderplot, aes(x = Group, y= Percentage))+
  facet_grid(. ~ Study)+
  scale_fill_manual(values = c("#499483", "#444c5c"))+
  #scale_fill_manual(values = c("#005073","#428bca"))+
  scale_alpha_discrete(range = c(0.5,1), 
                       guide = "none")+
  geom_col(aes(fill = Gender))+
  #alpha = Samplestatistic))+
  labs(x= "",
       y= "Percentage",
       alpha = "SampleStats")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        #strip.text = element_text(angle= 60)
        #legend.position = "none",
        text = element_text(size = 15),
        strip.text.x = element_blank())

plot.Gender


##################### Combine Age and Gender plots ##################### 

plot.GenderAge= plot_grid(plot.Age, plot.Gender, 
                          align="v",
                          axis = "rl", 
                          ncol =1, 
                          labels= "AUTO",
                          rel_heights = c(1.4,1))



####################### print PDFs
# Save Figures as PDF

pdf("data/figures/sup/SupplementalFigure1.pdf",
    width = 14,
    height = 6.5)
plot.GenderAge

dev.off()

#######################
library(WriteXLS)
  

supp = sample %>% 
  select(study)%>%
  rename(Study= study) %>% 
  left_join(ages %>% 
              rownames_to_column("Study") %>% 
              inner_join(gender %>%
                           rownames_to_column("Study"), by = "Study"
                         )
            )


WriteXLS(x= supp, 
         ExcelFileName = "data/paper_sup/gender_age.xlsx")
  