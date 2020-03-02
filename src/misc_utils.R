# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

library(tidyverse)

#' This function plots a fast MDS
#' 
#' @param FuncMatrix gene expression matrix with genes in rows and samples in columns
#' @param Targets data frame with 2 columns: Sample and HeartFailure, the order of Sample
#' has to be the same as the column order in ExpMat
#' @return an MDS plot
myMDS = function(FuncMatrix, Targets, main = "Title"){ 
  #So that samples are in the same order in the columns of the matrix
  library(limma)
  library(ggplot2)
  #and the id column of Targets
  
  MDS_res = plotMDS(FuncMatrix, plot = FALSE)
  
  myMDS_DF = data.frame(dim1 = MDS_res$x[Targets$Sample],
                        dim2 = MDS_res$y[Targets$Sample],
                        HeartFailure = Targets$HeartFailure)
  
  cbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  MDS_plot = ggplot(myMDS_DF,aes(x = dim1, y = dim2)) + geom_point(aes(color = factor(HeartFailure))) +
    theme_bw() + scale_colour_manual(values=cbPalette) + theme(panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank()) +
    ggtitle(main)
  
  print(MDS_plot)
}


#' This function calculates moderated t-stats using limma for experiments
#' in the metaheart object
#' 
#' @param ExpMat gene expression matrix with genes in rows and samples in columns
#' @param Targets data frame with 2 columns: Sample and HeartFailure, the order of Sample
#' has to be the same as the column order in ExpMat
#' @return a tibble with limma stats

run_HFlimma = function(ExpMat, Targets){
  library(limma)
  
  #Adjust a linear model to the expression data
  f = factor(Targets$HeartFailure, levels= c("yes","no"))
  design = model.matrix(~0+f)
  colnames(design) = c("yes","no")
  fit = lmFit(ExpMat, design)
  
  #Define contrasts
  cont.matrix = makeContrasts(HF_sign = yes-no,
                               levels=design)
  
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  
  #obtain differentially expressed genes
  HF_results = as.data.frame(topTable(fit2,adjust.method = "BH",number = Inf)) %>% 
    tibble::rownames_to_column() %>%
    arrange(desc(abs(t))) %>% as_tibble()
  
  colnames(HF_results)[1] = "ID"
  
  return(HF_results)
}

#' This function summarizes the results of msviper in an useful tibble
#' 
#' @param msviper_res results from an msviper run
#' @return a tibble with ordered viper stats


msviper_summary = function(msviper_res){
  msviper_tibble = as_tibble(data.frame(RegulonName = names(msviper_res$es$nes),
                                        NES = msviper_res$es$nes,
                                        Size = msviper_res$es$size[names(msviper_res$es$nes)],
                                        pvalue = msviper_res$es$p.value[names(msviper_res$es$nes)],
                                        adj_pvalue = p.adjust(msviper_res$es$p.value[names(msviper_res$es$nes)],
                                                              "BH"),
                                        stringsAsFactors = F)) %>% arrange(desc(abs(NES)))
  
}


#' This function calculates TF activities from a contrast in a METAheart study
#' 
#' @param x a study inside the METAheart object
#' @return a tibble with ordered viper stats


get_TFactivities = function(x){
  library(viper)
  load("HGEX_data/funcomics/BEST_viperRegulon.rdata") #dorothea
  
  dorothea_regulons = viper_regulon[-grep("_[ED]",
                                          names(viper_regulon))]
  
  dorothea_rank = x$HF_limma$t
  names(dorothea_rank) = x$HF_limma$ID
  dorothea_rank = sort(dorothea_rank)
  
  names(dorothea_regulons) = unlist(lapply((strsplit(names(dorothea_regulons),"_")), 
                                           function(x) x[1]))
  
  dorothea_results = msviper_summary(msviper(dorothea_rank,
                                             dorothea_regulons,
                                             minsize = 5,
                                             verbose = FALSE))
  
  return(dorothea_results)
  
}


#' This function performs GSEA to a contrast of a METAheart study
#' 
#' @param x a study inside the METAheart object
#' @return a tibble with ordered GSEA stats


get_GSEA = function(x){
  library(fgsea)
  load("HGEX_data/funcomics/TEDDY_geneSets.ro") #Gene set collection
  
  red_genesets = TEDDY_geneSets[c("MSIGDB_CANONICAL","MSIGDB_HMARKS",
                                 "MSIGDB_GO_CELLCOMP","MSIGDB_GO_BIOLPROC",
                                 "MSIGDB_GO_MOLFUNC","MSIGDB_CHEM-GEN")]
  
  gsea_rank = x$HF_limma$t
  names(gsea_rank) = x$HF_limma$ID
  #gsea_rank = sort(gsea_rank)
  
  set.seed(1234) # fgsea unstable results
  
  GSEA_results_complete = lapply(red_genesets, function(y){
    GSEA_results = fgsea(pathways = y, stats = gsea_rank,
                         minSize = 15, maxSize = 500, 
                         nperm = 1000) %>% as_tibble() %>% 
      arrange(desc(abs(NES)))
    
  })
  
  GSEA_results_complete = enframe(GSEA_results_complete) %>% 
                          unnest() %>% 
                          arrange(NES)
  
  return(GSEA_results_complete)
  
}

#' This function calculates PROGENy scores from a contrast in a METAheart study
#' and significance via a permutation approach
#' 
#' @param x a study inside the METAheart object
#' @return a tibble with ordered PROGENy stats


get_PROGENy_scores = function(x){
  library(progeny)
  model = as.matrix(read.csv(file = "HGEX_data/funcomics/ExtendedProgeny.txt",
                             row.names = 1)) #Extended PROGENy
  
  progeny_rank = x$HF_limma$t
  names(progeny_rank) = x$HF_limma$ID
  progeny_rank = sort(progeny_rank)
  
  progeny_rank_mat = as.matrix(progeny_rank)
  colnames(progeny_rank_mat) = "study"
  
  meta_progeny = t(progeny(progeny_rank_mat))
  
  #Significance: permutation
  permutation_prog = sapply(1:1000, function(x){
    rnd_order = progeny_rank_mat
    rownames(rnd_order) = sample(rownames(rnd_order))
    progeny(rnd_order)
  })
  
  rownames(permutation_prog) = rownames(meta_progeny)
  
  progeny_pvals = meta_progeny[,1]
  
  for(p in names(meta_progeny[,1])){
    
    metascore = meta_progeny[p,1]
    
    if(sign(metascore)>0){
      
      pval = sum(permutation_prog >= metascore)/length(permutation_prog)
      
    }else{
      
      pval = sum(permutation_prog <= metascore)/length(permutation_prog)
      
    }
    
    progeny_pvals[p] = pval
    
  }
  
  prog_res = tibble(progeny_scores = meta_progeny[,1], 
                    progeny_pvals, pathway = names(progeny_pvals)) %>%
    arrange(progeny_pvals)
  
  return(prog_res)
  
}

#' This function performs the meta analysis from a metaheart object
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param n_missing number of MAX NUMBER of studies that can miss a given gene
#' (inclusion criteria)
#' @return a sorted named vector of p-values from the Fisher Combined Test


run_fisher_meta = function(meta_list, n_missing = 3){
  library(survcomp)
  # Getting p-values from limma
  limma_pvals = get_all_limma(meta_list = meta_list, "adj.P.Val")
  
  # Use only genes that are present in all experiments (missing in n at most)
  limma_results_mat = limma_pvals[rowSums(is.na(limma_pvals))<=n_missing,]
  
  # Fisher combined test
  fisher_pvals = apply(limma_results_mat, 1, function(x){ 
    survcomp::combine.test(x, "fisher", na.rm = T)
  })
  
  fisher_pvals_adj = sort(p.adjust(fisher_pvals,"BH"))
  
  return(fisher_pvals_adj)
}

#' This function calculates a disease score matrix from an
#' expression matrix and a coefficient matrix  
#' 
#' 
#' @param CoefMatrix a coefficient matrix used to calculate the disease score:
#' genes in columns, predictor experiment name in columns
#' @param ExpressionMatrix group of genes to be used to calculate the disease score
#' @return disease score matrix, samples of index in rows, 
#'     predictor experiments in columns


getRiskIndex_matrix = function(CoefMatrix, ExpressionMatrix){
  
  RiskMatrix = apply(CoefMatrix, 2, function(t_vector){
    
    t_red = na.omit(t_vector)
    t_red = t_red[names(t_red) %in% rownames(ExpressionMatrix)]
    ExpressionMatrix_red = ExpressionMatrix[names(t_red),]
    RiskE = t_red %*% ExpressionMatrix_red
    
  } )
  
  RiskMatrix = scale(RiskMatrix)
  rownames(RiskMatrix) = colnames(ExpressionMatrix)
  
  
  return(RiskMatrix)
}


#' This function calculates a disease score and uses it 
#' to calculate pairwise classification performance
#' 
#' @param Experiment_List dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param genes group of genes to be used to calculate the disease score
#' @param limma_t_mat a coefficient matrix used to calculate the disease score:
#' genes in columns, predictor experiment name in columns
#' @return a named list, with predicted experiment as an index: 
#' RiskMatrix contains the disease score matrix, samples of index in rows, 
#'     predictor experiments in columns
#' SingleAUC contains a vector of AUCs coming of predicting index with a 
#'     single disease score (one column of RiskMatrix)
#' AUC_All contains a single AUC coming of predicting index with a consensus
#'     disease score. (not including index study if it appears in the RiskMatrix)


getRisk_Stats_v2 = function(Experiment_List,genes,limma_t_mat){
  library(ROCR)
  Experiments = names(Experiment_List)
  names(Experiments) = Experiments
  Risk_Stats = lapply(Experiments, function(E){
    
    E_GEX = Experiment_List[[E]]$GEX
    RiskMatrix =  getRiskIndex_matrix(CoefMatrix = limma_t_mat[genes,],
                                      ExpressionMatrix = E_GEX) # Gets risk matrix
    
    SingleAUC = apply(RiskMatrix,2, function(Risk_vect){ #Each column is one estimation of the risk index
      Risk_vect = Risk_vect[Experiment_List[[E]]$TARGETS$Sample]
      RiskDF = data.frame("RiskIX" = Risk_vect,
                          "HeartFailure" = ifelse(Experiment_List[[E]]$TARGETS$HeartFailure == "yes",1,0),
                          stringsAsFactors = F)
      
      ROC_E = prediction(RiskDF$RiskIX, RiskDF$HeartFailure) #Evaluate classification
      AUC_E = performance(ROC_E, measure = "auc")@y.values[[1]]
      return(AUC_E)
    })
    
    # Consensus score is calculated by taking the mean of each disease scores, without considering
    # the t values of the evaluated experiment
    ALL_Risk = rowMeans(RiskMatrix[Experiment_List[[E]]$TARGETS$Sample,colnames(RiskMatrix)!=E])
    
    RiskDF = data.frame("RiskIX" = ALL_Risk,
                        "HeartFailure" = ifelse(Experiment_List[[E]]$TARGETS$HeartFailure == "yes",
                                                1,0),
                        stringsAsFactors = F)
    
    ROC_All = prediction(RiskDF$RiskIX, RiskDF$HeartFailure)
    AUC_All = performance(ROC_All, measure = "auc")@y.values[[1]]
    
    RiskIX_results = list("RiskMatrix" = RiskMatrix,
                          "SingleAUC" = SingleAUC, #AUCs of predicting index experiment with experiments in t_mat
                          "AUC_All" = AUC_All) #AUCs of predicting index experiment with consensus of t_mat
    
    return(RiskIX_results)
    
  })
  return(Risk_Stats) #Index identical to Experiment_List
}

#' This function calculates a disease score and uses it 
#' to calculate pairwise classification performance
#' 
#' @param Experiment_List dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param genes group of genes to be used to calculate the disease score
#' @param limma_t_mat a coefficient matrix used to calculate the disease score:
#' genes in columns, predictor experiment name in columns
#' @return a named list, with predicted experiment as an index: 
#' RiskMatrix contains the disease score matrix, samples of index in rows, 
#'     predictor experiments in columns
#' SingleAUC contains a vector of AUCs coming of predicting index with a 
#'     single disease score (one column of RiskMatrix)
#' AUC_All contains a single AUC coming of predicting index with a consensus
#'     disease score. (not including index study if it appears in the RiskMatrix)


pairwise_ds = function(experiments, meta_list, ngenes =100, t_matrix){
  
  bind_rows(lapply(experiments, function(e){ #means that we get list of lists
    e_genes = (meta_list[[e]]$HF_limma[[1]])[1:ngenes] #Defining index-specific signature (predictor)
    e_genes = e_genes[e_genes %in% rownames(t_matrix)]
    e_ds_results = getRisk_Stats_v2(Experiment_List = meta_list,
                                    limma_t_mat = t_matrix,
                                    genes = e_genes) #Getting risk Stats
    
    AUC_n_louAUC = bind_rows(lapply(e_ds_results, function(Expment){
      #For each experiment in list (predicted) retrieve
      #performance of predictor Experiment
      return(list("single"= Expment$SingleAUC[e], #How predictor('s t values) predicted each experiment
                  "lou" = Expment$AUC_All)) #How predictor genes made the consensus better (all t's)
    }),.id = "PredictedExperiment")
    
    return(AUC_n_louAUC)
  }),.id = "PredictorExperiment") %>% arrange(PredictorExperiment,PredictedExperiment)
}

#' This function calculates the enrichment score of the up and downregulated genes
#' of the top n genes of each study, in the gene-stats of all studies.
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param ngenes number of top genes to be used
#' @return a unified tibble with mean abs(ES) per pairwise comparisons


pairwise_ES = function(meta_list, ngenes = 200){
  library(fgsea)
  
  study_deg_list_up = lapply(meta_list, function(x){
    deg = dplyr::slice(x$HF_limma,
                       1:ngenes) %>%
      filter(t >0)
    return(deg[[1]])
  })
  
  study_deg_list_down = lapply(meta_list, function(x){
    deg = dplyr::slice(x$HF_limma,
                       1:ngenes) %>%
      filter(t <0)
    return(deg[[1]])
  })
  
  up_ES = lapply(experiments, function(x){
    
    stat_rank = meta_list[[x]][["HF_limma"]][["t"]]
    names(stat_rank) = meta_list[[x]][["HF_limma"]][["ID"]]
    stat_rank = sort(stat_rank)
    
    up_row = as_tibble(fgsea(pathways = study_deg_list_up,
                             stats = stat_rank,nperm = 1000)) %>%
      dplyr::select(pathway,ES)
  })
  
  up_ES = up_ES %>% 
    enframe("Reference") %>% unnest()
  
  colnames(up_ES) = c("Reference","DEG","ES")
  
  up_ES = mutate(up_ES, direction = "up")
  
  down_ES = lapply(experiments, function(x){
    
    stat_rank = meta_list[[x]][["HF_limma"]][["t"]]
    names(stat_rank) = meta_list[[x]][["HF_limma"]][["ID"]]
    stat_rank = sort(stat_rank)
    
    up_row = as_tibble(fgsea(pathways = study_deg_list_down,
                             stats = stat_rank,nperm = 1000)) %>%
      dplyr::select(pathway,ES)
  })
  
  down_ES = down_ES %>% 
    enframe("Reference") %>% unnest()
  
  colnames(down_ES) = c("Reference","DEG","ES")
  
  down_ES = mutate(down_ES, direction = "down")
  
  all_ES = bind_rows(down_ES, up_ES)
  
  all_ES = mutate(all_ES, number_genes = ngenes)
  
  return(all_ES)
}

#' This function calculates an ANOVA of each feature in a PCA matrix and a covariate in target file
#' 
#' @param numeric_matrix samples in columns, features in rows
#' @param targets a tibble that contains at least the following columns: Sample, factor_a_name
#' @param factor_a string with the name of the factor to associate with each feature
#' @return a data frame with all anova_stats results for each feature in numeric_matrix

run_anovastats_single = function(numeric_matrix, targets, 
                                 factor_a = "ExpID",
                                 pval = 0.05){
  library(sjstats)

  pca_anova_study = apply(numeric_matrix, 1, function(x, targets){
    pc_i = x
    factor_a_vect = factor(targets[[factor_a]])
    gene_aov = aov(pc_i ~ factor_a_vect)
    aov_stats = sjstats::anova_stats(gene_aov) 
    
  },targets = targets) %>% 
    bind_rows(.id = "PC") %>% 
    as_tibble() %>% 
    group_by(PC) %>% 
    gather(stats,value,-(PC:term)) %>% 
    spread(term,value) %>%
    ungroup()
  
  # here we identify PCs associated with the study label
  pcs_study = pca_anova_study %>% 
    filter(stats == "p.value" & factor_a_vect < pval)
  
  return(pcs_study)
}


#' This function calculates a 2-way ANOVA with interactions
#' 
#' @param numeric_matrix samples in columns, features in rows
#' @param targets a data.frame that contains at least the following columns: Sample, factor_a_name, factor_b_name
#' @param factor_a string with the name of the first factor
#' @param factor_b string with the name of the second factor
#' @return a data frame with all anova_stats results for each feature in numeric_matrix


run_anovastats = function(numeric_matrix, targets, 
                          factor_a = "HeartFailure", factor_b = "ExpID"){
  library(sjstats)
  library(fastmatch)
  
  anova_res = apply(numeric_matrix, 1, function(x, targets){
    
    gene_i = x
    factor_a_vect = factor(targets[[factor_a]])
    factor_b_vect = factor(targets[[factor_b]])
    gene_aov = aov(gene_i ~ factor_a_vect * factor_b_vect)
    aov_stats = sjstats::anova_stats(gene_aov) 
    
  },targets = targets)  %>% bind_rows(.id = "gene") %>% as_tibble() %>% 
    group_by(gene) %>% gather(stats,value,-(gene:term))  %>% spread(term,value) %>%
    ungroup()
  
  colnames(anova_res)[fmatch(c("factor_a_vect","factor_b_vect","factor_a_vect:factor_b_vect"),
                             colnames(anova_res))] = c(factor_a,factor_b,
                                                       paste(factor_a,factor_b,sep="_"))
  
  return(anova_res)
  
}