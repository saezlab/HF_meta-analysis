# MIT License

# Copyright (c) [2020] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description:  Analysis and plotting of the Human Protein Atlas data to select biomarker candidates


# Load data and libraries ---------------------------------------------------------------------

library(cowplot)
library(tidyverse)

# HPA 
protatlas =as_tibble(read.delim("data/paper_sup/normal_tissue.tsv"))

# the HF-CS:
meta_rank = read.csv(file =  "data/shiny/meta_analysis_summary.txt",sep = "\t",
                    header = T, stringsAsFactors = F) %>%
  mutate(log10pval = -log10(fisher_pvalue)) %>%
  arrange(desc(log10pval)) %>%
  mutate(rank = seq(1, dim(.)[1]))

# fgsea results from validation_plotting.R 
fgseaRes= readRDS(file = "data/figure_objects/Validation_GSEA_results.rds")  

#biomarker candidates form our analysis:
prot_manifest =readRDS("data/figure_objects/proteomics_manifest.rds") %>% 
  inner_join(meta_rank %>% rename(hgnc_symbol = gene), by = "hgnc_symbol") %>%
  mutate(leading = hgnc_symbol %in% fgseaRes$leadingEdge[[1]],
         quadrant= sign(logor) == sign(trnscrpt_evid_t),
         colored = (quadrant==T) & (leading == T),
         metaTop500= hgnc_symbol %in% meta_rank$gene[1:500] )

# Analysis HPA - SF panel A -------------------------------------------------------------------------------

#merge data HPA and biomarker data:
candidates = protatlas %>% 
  filter(Tissue == "heart muscle") %>% 
  rename(hgnc_symbol = Gene.name) %>% 
  right_join(prot_manifest) %>%
  left_join(meta_rank)


#filter for interesting candidates (as ploted in Fig5 in the main manuscript)
plot_cand = candidates %>% 
  filter(colored == T) %>% # colored indicates leading edge genes that agree with the direction of regulation
  mutate(Level = factor(Level, levels =c("High", "Medium","Low",  "Not detected")),
         Reliability = factor(Reliability, levels =c("Enhanced", "Supported","Approved",  "Uncertain")))%>%
  arrange(rank)

#filter for candidates that are detected by the HPA
cand_list =plot_cand %>% 
  filter(Reliability != "Uncertain", Level != "Not detected") %>% 
  arrange(rank) %>% 
  pull(hgnc_symbol)

#create a colour vector for plot_cand that are in cand_list
con <- ifelse(plot_cand %>% 
                arrange(rank) %>%
                drop_na() %>%
                pull(hgnc_symbol) %in% cand_list, 'darkblue', 'black')


# plot protein expression in heart tissue- SF panel A 
plot.prot.HPA= ggplot(plot_cand %>% drop_na, aes(x= reorder(hgnc_symbol,rank), 
                                                 size= Level, 
                                                 color= Reliability,
                                                 y= rank ))+
  geom_count()+
  scale_color_manual(values = c("#660000","#CC0000","#FF3333","#FF9999", "#FFE5CC"))+
  scale_size_manual(values = c(6,4.5,3,1))+
  theme_minimal()+
  labs(x= "Biomarker candidates", 
       y= "HF-CS rank")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ggtitle("Protein expression in heart muscle by HPA ")+
  scale_y_reverse()+
  theme(axis.text.x = element_text(color = con))


# HPA tissue specificity - SF panel B ---------------------------------------------------------

# count how many times each candidate was reported in other tissues than heart
# HPA data is select for non cardiac tissue and non "uncertain" information
# as multiple cell types are reported per tissue, the highest level of the ordinal variable "Level"
# is selected for a tissue.

tissue_count= sapply(cand_list, function(x){
  y = protatlas %>% filter(Gene.name ==x) %>% 
    filter(#Level != "Not detected",
      #Level != "Low",
      Reliability != "Uncertain", 
      Tissue != "heart muscle") %>%
    group_by(Tissue) %>%
    mutate(Level= factor(Level, levels = c("Not detected", "Low", "Medium","High"), ordered = T))%>%
    mutate(high_Level = max(Level)) %>%
    ungroup()%>%
    distinct(Gene.name, Tissue, high_Level)
  
  c("ntissue"= length(unique(y$Tissue)),table(y$high_Level))

  })

# Score is now calculated by assesing the proportion of tissue number with High and Medium expression levels 
# of all measured tissues. 

tissue_spec = t(tissue_count) %>% 
  as_tibble() %>% 
  mutate("Biomarker" = rownames(t(tissue_count))) %>%
  mutate(score= (Low + `Not detected`) / ntissue)

plot_tissue = ggplot(tissue_spec, aes(x= reorder(x= Biomarker,-score), y= score))+
  geom_point(size =3)+
  theme_minimal()+
  labs(x = "Biomarker candidates expressed in heart muscle",
       y= "Cardiac tissue specificity")+
  #scale_y_reverse()+
  ggtitle("Tissue specificity based on the HPA")+
  theme(axis.text.x = element_text(color = "darkblue"))
plot_tissue




# Combine Plots -------------------------------------------------------------------------------

### combin plots (supp figure X)

plot_hpa = plot_grid(plot.prot.HPA, plot_tissue, nrow= 2, labels = "AUTO", 
                     rel_heights = c(1,0.8))

pdf("data/figures/sup/SupplementalFigure_HPA.pdf",
    width = 7,
    height = 6)
plot_hpa
dev.off()

