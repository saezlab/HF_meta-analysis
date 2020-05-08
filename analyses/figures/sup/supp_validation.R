# MIT License

# Copyright (c) [2020] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: Plotting supplement figure, displayign additional validation results
# 1. fetal transcriptome results
# 2. development hf biomarker

library(tidyverse)
library(fgsea)
library(cowplot)
library(patchwork) 
library(ggrepel)

fgseaRes= readRDS("data/figure_objects/Validation_GSEA_results.rds")

meta_rank = read.csv(file =  "data/shiny/meta_analysis_summary.txt",sep = "\t",
                     header = T, stringsAsFactors = F)

meta_rank = meta_rank %>%
  mutate(log10pval = -log10(fisher_pvalue)) %>%
  arrange(desc(log10pval)) %>%
  mutate(rank = seq(1, dim(.)[1]))




#################### Fetal plots (Genes and TFs) second study GSE52601
##GENE
fetalDEA= readRDS(file ="data/figure_objects/fetalDEgenes_GSE52601.rds") %>% 
  inner_join(meta_rank, by = "gene") %>% 
  filter(adj.P.Val< 0.05) %>%
  mutate(leading = gene %in% fgseaRes$leadingEdge[[3]],
         quadrant= sign(t) == sign(mean_t),
         colored = (quadrant== T) & (leading == T), 
         metaTop500 = (gene %in% meta_rank$gene[1:500])) %>% 
  arrange(desc(fisher_pvalue)) 

plot.fetal.gene.GSE52601 = ggplot(data= fetalDEA, aes(x= t, y= mean_t))+
  geom_point(aes(color =colored ))+
  scale_color_manual(values =c("#d0d0e1","#000000"))+
  #ggtitle("fetal transcriptome vs HF consensus signature")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("t-value (GSE52601)")+
  ylab("Mean t-value in HF")+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  geom_vline(xintercept= 0, color = "grey", linetype = 2)+
  geom_label_repel(aes(label=ifelse( (metaTop500 == T) &
                                       (-log10(adj.P.Val)> 4.3) & 
                                       (colored == T)
                                     ,gene
                                     ,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size =3)
plot.fetal.gene.GSE52601

##TF
plot.fetal.TFs.data= readRDS(file="data/figure_objects/fetalTFs_GSE52601.rds")

plot.fetal.TFs.GSE52601 = ggplot(data= plot.fetal.TFs.data, aes(x= NES.fetal,y = NES.meta))+
  geom_point(aes(color =colored ))+
  scale_color_manual(values =c("#d0d0e1","#000000"))+
  #ggtitle("Plasma proteomics vs HF consensus signature")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("TF activity in GSE52601 (NES)")+
  ylab("TF activity in HF (NES)")+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  geom_vline(xintercept= 0, color = "grey", linetype = 2)+
  geom_label_repel(aes(label=ifelse((metasig == T) &
                                      (colored == T)
                                    ,RegulonName,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size = 3)


# Development proteome
prot_dev = readRDS("data/figure_objects/proteomics_dev.rds")

prot_dev= prot_dev %>% inner_join(meta_rank %>% rename(hgnc_symbol = gene), by = "hgnc_symbol") %>%
  mutate(leading = hgnc_symbol %in% fgseaRes$leadingEdge[[2]],
         quadrant= sign(loghr) == sign(trnscrpt_evid_t),
         metaTop500= hgnc_symbol %in% meta_rank$gene[1:500],
         colored = (metaTop500==T) & (quadrant== T))


plot.proteomics.earlyHF = ggplot(data= prot_dev, aes(x= loghr,y = trnscrpt_evid_t))+
  geom_point(aes(color =colored))+
  scale_color_manual(values =c("#d0d0e1","#000000"))+
  #ggtitle("Plasma proteomics vs HF consensus signature")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Early HF plasma proteome log2(HR)")+
  ylab("Mean t-value in HF")+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  geom_vline(xintercept= 0, color = "grey", linetype = 2)+
  geom_label_repel(aes(label=ifelse( (metaTop500 == T) &
                                       (colored == T)
                                     ,hgnc_symbol,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size =3)

# add disease score classifier for fetal samples (output from ds_fetal_external_studies.R)
#ds_fetal = readRDS("data/figure_objects/ds_fetal_plot.rds")
ds_fetal = readRDS("data/figure_objects/ds_fetal_experiment.rds") %>%
  enframe() %>% unnest()

plot.ds.fet = ggplot(ds_fetal,
                     aes(x=name, y=Risk_Score, color=HeartFailure)) +
  geom_hline(yintercept = 0,
             color = "grey",
             linetype = "dashed")+
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values=c("005073","#499483"),
                      name = "",
                      labels = c("CT", "fetal") ) + 
  labs(y= "Disease score",
       x= "",
       color = "Fetal Sample")+
  theme_minimal()+
  ylim(c(-2.5,2.1))+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_line(size =0.5),
        axis.text = element_text(size= 11.5),
        axis.title.y= element_text(size =13))


## combine plots

suppfigure= plot_grid(plot.ds.fet, 
                      plot.fetal.gene.GSE52601, 
                      plot.fetal.TFs.GSE52601,
                      plot.proteomics.earlyHF,
                      ncol = 2,
                      rel_widths = c(1,1.2),
                      labels= "AUTO")

pdf("data/figures/sup/SupplementalFigure13.pdf",
    width = 9,
    height = 9)
suppfigure
dev.off()

