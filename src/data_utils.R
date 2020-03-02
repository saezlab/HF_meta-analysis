# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# General functions used in all analysis

library(tidyverse)

#' This function combines a group of tibbles that contain the same information
#' but come from different experiments
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment _> list of tables of information
#' @param index_name name of the tables to combine uniformly
#' @return joint tibble with an extra variable containing the name of the experiment


get_tibble_union = function(meta_list, index_name){
  experiments = names(meta_list)
  meta_tible = tibble()
  
  for(e in experiments){
    meta_list[[e]][[index_name]] = mutate(meta_list[[e]][[index_name]],
                                          ExpID = e)
    
    meta_tible = dplyr::bind_rows(meta_list[[e]][[index_name]],
                                  meta_tible)
  }
  
  return(meta_tible)
}

#' This function applies the tibble union to get differential expression analysis results of all studies in the
#' meta-analysis
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param limma_column define the statistic from limma results to be fetched
#' @return a matrix with genes in the rows, experiments in the columns and statistics in the cells


get_all_limma = function(meta_list, limma_column){
  
  sel_cols =  c("ID","ExpID", limma_column)
  
  limma_results = get_tibble_union(meta_list,"HF_limma") %>% 
    dplyr::select(sel_cols) %>% 
    spread(sel_cols[2],sel_cols[3])
  
  limma_results_mat =  as.matrix(limma_results[,-1])
  rownames(limma_results_mat) = limma_results[[1]]
  
  return(limma_results_mat)
}

#' This function applies the tibble union to get differential expression analysis results of all studies in the
#' meta-analysis
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param id_column define the id column from the slot selected to be fetched
#' @param stat_column define the statistic column from the slot selected to be fetched
#' @param slot which results to merge
#' @return a matrix with genes in the rows, experiments in the columns and statistics in the cells
#' 
get_all_stats = function(meta_list, id_column, stat_column, slot){
  
  sel_cols =  c(id_column,"ExpID", stat_column)
  
  slot_results = get_tibble_union(meta_list,slot) %>% 
    dplyr::select(sel_cols) %>% 
    spread(sel_cols[2],sel_cols[3])
  
  slot_results_mat =  as.matrix(slot_results[,-1])
  rownames(slot_results_mat) = slot_results[[1]]
  
  return(slot_results_mat)
}

#' This function generates a unified gene expression matrix of all experiments defined in a
#' target data frame/tibble
#' 
#' @param meta_list dictionary that contains all the information of studies. 
#'     Usually METAheart object: Experiment -> list of tables of information
#' @param gex_key where in the meta_list is the matrix you want to merge
#' @param complete_targets unified target file of all experiments that will be merged
#' @return a matrix with genes in the rows, general ids in the columns 
#' (exp_id + original sample name) and expression in the cells


# Get unified gene expression matrix
get_complete_gex = function(meta_list, gex_key = "GEX",
                            complete_targets){
  #METAheart object is a 
  
  all_genes = unique(unlist(lapply(meta_list,
                                   function(x) rownames(x[[gex_key]]))))
  
  gex_union = matrix(NA,nrow = length(all_genes),ncol = nrow(complete_targets))
  
  colnames(gex_union) = complete_targets$grl_id
  rownames(gex_union) = all_genes
  
  for(e in unique(complete_targets$ExpID)){
    fdf = filter(complete_targets,ExpID == e)
    gex = (meta_list[[e]])[[gex_key]]
    genes = rownames(gex)[rownames(gex) %in% all_genes]
    gex = gex[genes,fdf$Sample]
    gex_union[genes,fdf$grl_id] = gex
  }
  
  return(gex_union)
  
}

## Function to group Dorothea regulons. 
## Input: A data frame containing Dorothea regulons, as stored in 
## https://github.com/saezlab/ConservedFootprints/tree/master/data
## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}
