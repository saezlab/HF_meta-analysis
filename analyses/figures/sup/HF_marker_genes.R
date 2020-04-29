# Author: Jan Lanzer
# Date: Nov 2019
# Checking Geneexpression of important HF genes per study

library(tidyverse)

METAheart= readRDS("HGEX_data/METAheart.rds")
Experiments = names(METAheart)

### 1) Define the set of genes that should be up or down regulated in HF

HFgenes_up = c("NPPA","NPPB","MYH7","MME", "COL1A1","COL3A1", "POSTN", "MMP2", "MMP9")
HFgenes_dn = c("MYH6", "ATP2A2","SLC2A1", "KCNH2", "TNNT2")

### 2) Analyzing the METAheart study for those genes
# creating the dataframe to store the t values of the defined gene lists

HFgenes = data.frame(matrix(NA, ncol = length(Experiments), nrow = length(c(HFgenes_up, HFgenes_dn))))
colnames(HFgenes) = Experiments
rownames(HFgenes) = c(HFgenes_up, HFgenes_dn)

# loop over studies and genes to copy t-values into the dataframe
for (study in names(METAheart)){
    studygenes = METAheart[[study]]$HF_limma %>%
      filter(ID %in% c(HFgenes_dn, HFgenes_up))
  for (gene in studygenes$ID){
    HFgenes[gene,study] = studygenes %>% filter(ID == gene) %>% 
      select(t)
    }
  }

# tidy formatting
HFgenes_tidy = HFgenes %>% 
  rownames_to_column("gene") %>% 
  as.tibble() %>% 
  gather(Experiments, key = Study, value = "t")%>% 
  mutate(gene = factor(gene, levels = c(HFgenes_up, HFgenes_dn)))

## plot results
plot.HF = ggplot(data = HFgenes_tidy, aes(x= gene, y= t))+
  geom_boxplot()+
  geom_point(aes(color =Study), 
             size = 4, 
             alpha = 0.6)+
  theme_classic()+
  labs(x="HF marker genes",
       y = "t-value")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=12))+
  geom_hline(yintercept = 0, 
             color = "grey",
             linetype = 2)+
  geom_vline(xintercept = (length(HFgenes_up)+ 0.5),
             color = "black",
             linetype =1)+
  ggtitle("HF-marker gene expression (t-values) for all studies in metaheart project")


print(plot.HF)


### 3) Save figure

pdf("figures_pdfs/suppFig4_HF_markergenes.pdf",
    width = 12,
    height = 5)

plot.HF

dev.off()

