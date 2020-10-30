## what are genes undetected by single studies?
library(tidyverse)
library(cowplot)
MetaHeart= readRDS(file = "data/METAheart.rds")
HF_CS= read_delim("data/shiny/meta_analysis_summary.txt", delim= '\t')
#source("src/data_utils.R") #general functions 
#source("src/misc_utils.R")

#short= MetaHeart[1:2]

#function to retrieve all DE genes based on alpha

get_single_DE = function(Meta, alpha){
  lapply(Meta, function(x) {
    x$HF_limma %>% 
      filter(adj.P.Val < alpha) %>% 
      pull(ID)
  })
}


# set alpha
alpha= 0.1

# gent single study DE genes as a list
DE= get_single_DE(MetaHeart, alpha)

# plot how many DE genes per study: 
counts= sapply(DE, length)
counts= tibble(Study= names(counts), DE_genes= counts)

p1= ggplot(counts, aes(x= reorder(Study, DE_genes), y= DE_genes))+
  geom_col()+
  ggtitle(paste("DE genes at alpha <", alpha))+
  theme_minimal()+
  theme(text = element_text(size = 12), 
        axis.text.x = element_text(angle= 30))+
  labs(x= "Studies")

# count for every gene in HF_CS how many studies reported that gene 
study_counts= enframe(table(unlist(DE)), name= "gene", value= "studies_reported")

HF_CS= HF_CS %>% 
  left_join(study_counts, by= "gene") %>%
  mutate(studies_reported= as.numeric(studies_reported)) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(rank = c(1:dim(HF_CS)[1])) #%>% 
  #mutate(studies_reported= as.integer(studies_reported))
  
range(HF_CS$studies_reported)

p2 = ggplot(HF_CS, aes(x= studies_reported))+
  geom_histogram()+
  theme_minimal()


p3 = ggplot(HF_CS, aes(x= rank, y= studies_reported))+
  geom_point(size = 0.5)+
  theme_minimal()+
  scale_y_continuous(breaks = c(0:10))



HF_CS 

plot_grid(p1, p2, p3, 
          nrow= 2, 
          labels = "AUTO")


p4 = ggplot(HF_CS[1:500,], aes(x= studies_reported))+
  geom_histogram()+
  theme_minimal()+
  scale_x_continuous(breaks= c(0:10), limits = c(0,10))+
  labs(x= "number of studies that reported a gene at adj. p-value < 0.1",
       y= "gene count")
p4
 ggplot(HF_CS[501:1500,], aes(x= studies_reported))+
  geom_histogram()+
  theme_minimal()+
  scale_x_continuous(breaks= c(0:10), limits = c(0,10))+
  labs(x= "number of studies that reported a gene at adj. p-value < 0.1",
       y= "gene count")


p5 = ggplot(HF_CS[1:500,], aes(x= rank, y= studies_reported))+
  geom_point(size = 0.5)+
  theme_minimal()+
  scale_y_continuous(breaks = c(0:10))


plot_grid(p4, p5,  
          nrow= 1, 
          labels = "AUTO")


HF_CS[1:500,] %>% arrange(studies_reported, fisher_pvalue)

#### plot top 500 distributions: 

HF_CS2 = HF_CS %>% mutate(studies_reported= as.factor(studies_reported))
p6 = ggplot(HF_CS2[1:500,], aes(x= studies_reported))+
  geom_histogram(stat= "count")+
  scale_x_discrete(drop=F) +
  theme_minimal()+
  labs(x= "number of studies that reported a gene at adj. p-value < 0.1",
       y= "top 500 HF_CS")

p6
p7 = ggplot(HF_CS2[501:5000,], aes(x= studies_reported))+
  geom_histogram(stat= "count")+
  theme_minimal()+
  scale_x_discrete(drop=F) +
  #scale_x_continuous(breaks= c(0:10), limits = c(0,10))+
  labs(x= "Number of studies that reported a gene at adj. p-value < 0.1",
       y= "top 501-5000 HF_CS")

Plot_added_genes = plot_grid(p6+theme(axis.title.x = element_blank()) ,
                             p7,
                             nrow =2,
                             labels = "AUTO")

pdf(file = "data/figures/sup/DE_reported.pdf")
Plot_added_genes
dev.off()

unique(HF_CS$studies_reported)
gene= "TTC3"
studies= lapply(DE, function(x){
  gene %in% x
}) 
unlist(studies)

DE$Spurrell19
