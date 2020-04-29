  # MIT License
  
  # Copyright (c) [2020] [Jan D. Lanzer]
  # jan.lanzer@biquant.uni-heidelberg.de
  
  # Description: Plotting Validation approaches
  
  library(tidyverse)
  library(fgsea)
  library(cowplot)
  library(patchwork) 
  
  
  #################### 1. Perform Enrichment of all validation gene sets in HF consensus signature
  
  proteinsets= readRDS("HGEX_data/protein_genesets.rds")
  fetal_genesets_GSE = readRDS("HGEX_data/figure_objects/fealgenesetGSE52601.rds")
  fetal_gensets_P = readRDS("HGEX_data/figure_objects/fealgenesetPRJNA522417.rds")
  
  meta_rank = as_tibble(read.csv("HGEX_data/METArank_March2020.csv",
                                 header = T,
                                 sep= ",",
                                 stringsAsFactors = F) ) %>%
    mutate(log10pval = -log10(fisher_pvalue)) %>%
    arrange(desc(log10pval))
  
  gsea_rank_undir= meta_rank$log10pval
  names(gsea_rank_undir) = meta_rank$gene
  
  gsea_genesets = c(proteinsets,fetal_genesets_GSE,fetal_gensets_P)
  names(gsea_genesets) = c("plasma_prot_manifest_HF", 
                           "plasma_prot_early_HF", 
                           "fetal_GSE52601", 
                           "fetal_Spurrell19" )
  
  fgseaRes2 = fgsea(pathways=  gsea_genesets,
                         stats = gsea_rank_undir) %>% # update package, use higher permutations
    as_tibble()
  
  ## Plot enrichment pvals:
  
  plot.enrichment =ggplot(data= fgseaRes, aes(x= pathway,y= -log10(pval), fill= sign(NES)))+
    geom_bar(stat = "identity")+
    ylim(c(0,4))+
    coord_flip()+
    xlab("")+
    theme_minimal()+
    theme(text= element_text(size =12),
          legend.position = "none" )+
    geom_hline(yintercept= -log10(0.05), color = "grey", linetype = 2)
  
  
  #################### 2.Fetal plots (Genes and TFs)
  ##GENE
  fetalDEA= readRDS(file ="HGEX_data/figure_objects/fetalDEgenes_PRJNA522417.rds") %>% 
    inner_join(meta_rank, by = "gene") %>% 
    filter(adj.P.Val< 0.05) %>%
    mutate(leading = gene %in% fgseaRes$leadingEdge[[3]],
          quadrant= sign(t) == sign(mean_t),
          colored = (quadrant== T) & (leading == T), 
          metaTop500 = (gene %in% meta_rank$gene[1:500])) %>% 
    arrange(desc(fisher_pvalue)) 
  
  plot.fetal.gene = ggplot(data= fetalDEA, aes(x= t, y= mean_t))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("fetal transcriptome vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("t-value in fetal transcriptome")+
    ylab("Mean t-value in HF consensus")+
    geom_hline(yintercept = 0, color = "grey", linetype = 2)+
    geom_vline(xintercept= 0, color = "grey", linetype = 2)+
    geom_label_repel(aes(label=ifelse( (metaTop500 == T) &
                                         (-log10(adj.P.Val)> 8.9) & 
                                         (colored == T)
                                       ,gene
                                       ,"")),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     size =3)
  
  print(plot.fetal.gene)
  
  ##TF
  plot.fetal.TFs.data= readRDS(file="HGEX_data/figure_objects/fetalTFs_PRJNA522417.rds")
  
  plot.fetal.TFs = ggplot(data= plot.fetal.TFs.data, aes(x= NES.fetal,y = NES.meta))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("Plasma proteomics vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("TF activity in fetal samples (NES)")+
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
  
  plot.fetal.TFs
  
  #################### 3. Proteome plot
  
  # Manifest Proteome
  prot_manifest =readRDS("HGEX_data/proteomics_manifest.rds") %>% 
      inner_join(meta_rank %>% rename(hgnc_symbol = gene), by = "hgnc_symbol") %>%
    mutate(leading = hgnc_symbol %in% fgseaRes$leadingEdge[[1]],
           quadrant= sign(logor) == sign(trnscrpt_evid_t),
           colored = (quadrant==T) & (leading == T),
           metaTop500= hgnc_symbol %in% meta_rank$gene[1:500] )
  
  # grey + black dots: 1) significant hits in independent study
  # black dots: 2) part of leading edge , 3) same direction
    # labels: 4) in meta top 500
    
  plot.proteomics = ggplot(data= prot_manifest, aes(x= logor,y = trnscrpt_evid_t))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("Plasma proteomics vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("Proteomics log2(OR)")+
    ylab("Mean t-value in HF consensus")+
    geom_hline(yintercept = 0, color = "grey", linetype = 2)+
    geom_vline(xintercept= 0, color = "grey", linetype = 2)+
    geom_label_repel(aes(label=ifelse( (metaTop500 == T) &
                                         (colored == T)
                                       ,hgnc_symbol,"")),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     size =3)
  
  
  print(plot.proteomics)
  
  
  # Development proteome
  prot_dev = readRDS("HGEX_data/proteomics_dev.rds")
  dev_all_plt = ggplot(prot_dev, 
                       aes(x= loghr,
                           y = trnscrpt_evid_t,
                           label = hgnc_symbol)) +
    geom_point() + theme_minimal() + 
    geom_text_repel(size=3) + 
    theme(axis.title = element_text(size =12),
          axis.text = element_text(size =12),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black',
                                          size=1)) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    xlab("Proteomics log2(HR)") + ylab("Meta-analysis mean t-value")
  
  #### combine plots (fig6)
  
  plot4 = plot_grid(plot.proteomics, 
                    plot.fetal.TFs, 
                    ncol =1,
                    labels = c("C", "E"))
  blank_p <- plot_spacer() + theme_void()
  
  plot.enrichment2= plot_grid(blank_p, blank_p, plot.enrichment, ncol = 1, labels = c("A", "", "B"))
  plot1 = plot_grid(plot.enrichment2 ,plot.fetal.gene, ncol = 1, labels= c("","D"))
  
  plot5 = plot_grid( plot1, plot4, rel_widths = c(1.2,1))
  
  
  pdf("HGEX_figures/main/fig6.pdf",
      width = 10,
      height = 10)
  plot5
  dev.off()
  
  
  ########## create plots for supplemental figure
  
  
  #################### Fetal plots (Genes and TFs) second study GSE52601
  ##GENE
  fetalDEA= readRDS(file ="HGEX_data/figure_objects/fetalDEgenes_GSE52601.rds") %>% 
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
    xlab("t-value in fetal transcriptome")+
    ylab("Mean t-value in HF consensus")+
    geom_hline(yintercept = 0, color = "grey", linetype = 2)+
    geom_vline(xintercept= 0, color = "grey", linetype = 2)+
    geom_label_repel(aes(label=ifelse( (metaTop500 == T) &
                                         (-log10(adj.P.Val)> 4.5) & 
                                         (colored == T)
                                       ,gene
                                       ,"")),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     size =3)
  
  
  ##TF
  plot.fetal.TFs.data= readRDS(file="HGEX_data/figure_objects/fetalTFs_GSE52601.rds")
  
  plot.fetal.TFs.GSE52601 = ggplot(data= plot.fetal.TFs.data, aes(x= NES.fetal,y = NES.meta))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("Plasma proteomics vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("TF activity in fetal samples (NES)")+
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
  
  
  ## combine plots
  
suppfigure= plot_grid(plot.fetal.gene.GSE52601, plot.fetal.TFs.GSE52601,
                      ncol = 2,
                      labels= "AUTO",)
suppfigure

pdf("HGEX_figures/sup/fetal_comparison_GSE.pdf",
    width = 10,
    height = 5)
suppfigure
dev.off()

