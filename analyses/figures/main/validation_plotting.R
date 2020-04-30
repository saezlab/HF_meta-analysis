  # MIT License
  
  # Copyright (c) [2020] [Jan D. Lanzer]
  # jan.lanzer@biquant.uni-heidelberg.de
  
  # Description: Plotting Validation approaches
  # 1. Enrichment analysis of independent study results in HF consensus signature
  # 2. Fetal study Spurell19, A) DE genes, B) significant TFs
  # 3. Manifest HF plasma proteome, plot candidates
  
  library(tidyverse)
  library(fgsea)
  library(cowplot)
  library(patchwork) 
  library(ggrepel)
  
  
  #################### 1. Perform Enrichment of all validation gene sets in HF consensus signature
  
  proteinsets= readRDS("data/figure_objects/protein_genesets.rds")
  fetal_genesets_GSE = readRDS("data/figure_objects/fealgenesetGSE52601.rds")
  fetal_gensets_P = readRDS("data/figure_objects/fealgenesetPRJNA522417.rds")
  
  meta_rank = as_tibble(read.csv("data/METArank_March2020.csv",
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
  
  fgseaRes = fgsea(pathways=  gsea_genesets,
                         stats = gsea_rank_undir,
                   nperm= 1000,
                   minSize = 2,
                   gseaParam = 2) %>% 
              as_tibble()
  
  saveRDS(fgseaRes, file = "data/figure_objects/Validation_GSEA_results.rds")
              
  
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
  fetalDEA= readRDS(file ="data/figure_objects/fetalDEgenes_PRJNA522417.rds") %>% 
    inner_join(meta_rank, by = "gene") %>% 
    filter(adj.P.Val< 0.05) %>%
    mutate(leading = gene %in% fgseaRes$leadingEdge[[4]],
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
    xlab("t-value (Spurrell19)")+
    ylab("Mean t-value in HF")+
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
  plot.fetal.TFs.data= readRDS(file="data/figure_objects/fetalTFs_PRJNA522417.rds")
  
  plot.fetal.TFs = ggplot(data= plot.fetal.TFs.data, aes(x= NES.fetal,y = NES.meta))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("Plasma proteomics vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("TF activity (NES)(Spurrell19)")+
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
  prot_manifest =readRDS("data/figure_objects/proteomics_manifest.rds") %>% 
      inner_join(meta_rank %>% rename(hgnc_symbol = gene), by = "hgnc_symbol") %>%
    mutate(leading = hgnc_symbol %in% fgseaRes$leadingEdge[[1]],
           quadrant= sign(logor) == sign(trnscrpt_evid_t),
           colored = (quadrant==T) & (leading == T),
           metaTop500= hgnc_symbol %in% meta_rank$gene[1:500] )
  
  plot.proteomics = ggplot(data= prot_manifest, aes(x= logor,y = trnscrpt_evid_t))+
    geom_point(aes(color =colored ))+
    scale_color_manual(values =c("#d0d0e1","#000000"))+
    #ggtitle("Plasma proteomics vs HF consensus signature")+
    theme_minimal()+
    theme(legend.position = "none")+
    xlab("plasma_prot_manifest_HF log2(OR)")+
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
  # What is plotted? 
  # grey + black dots: significant hits in independent study
  # black dots: part of leading edge (from GSEA), same direction in consensus signature and 
  # and independent study
  # labels: hits that are rankend within top 500 genes in consensus signature

  #### combine plots (fig6)
  
  plot4 = plot_grid(plot.proteomics, 
                    plot.fetal.TFs, 
                    ncol =1,
                    labels = c("C", "E"))
  blank_p <- plot_spacer() + theme_void()
  
  plot.enrichment2= plot_grid(blank_p, blank_p, plot.enrichment, ncol = 1, labels = c("A", "", "B"))
  plot1 = plot_grid(plot.enrichment2 ,plot.fetal.gene, ncol = 1, labels= c("","D"))
  
  plot5 = plot_grid( plot1, plot4, rel_widths = c(1.2,1))
  
  pdf("data/figures/main/fig6.pdf",
      width = 10,
      height = 10)
  plot5
  dev.off()
  
  
